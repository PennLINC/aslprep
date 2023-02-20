"""BIDSLayout class."""
import os
import json
import re
from collections import defaultdict
from io import open
from functools import partial
from itertools import chain
import copy
import warnings
import sqlite3
import enum
from pathlib import Path

import sqlalchemy as sa
from sqlalchemy.orm import joinedload
from bids_validator import BIDSValidator

from ..utils import listify, natural_sort, make_bidsfile
from ..external import inflect
from .writing import build_path, write_contents_to_file
from .models import (Base, Config, BIDSFile, Entity, Tag)
from .index import BIDSLayoutIndexer
from .. import config as cf

try:
    from os.path import commonpath
except ImportError:
    def commonpath(paths):
        prefix = os.path.commonprefix(paths)
        if not os.path.isdir(prefix):
            prefix = os.path.dirname(prefix)
        return prefix

__all__ = ['BIDSLayout']


def parse_file_entities(filename, entities=None, config=None,
                        include_unmatched=False):
    """Parse the passed filename for entity/value pairs.

    Parameters
    ----------
    filename : str
        The filename to parse for entity values
    entities : list or None, optional
        An optional list of Entity instances to use in extraction.
        If passed, the config argument is ignored. Default is None.
    config : str or :obj:`bids.layout.models.Config` or list or None, optional
        One or more :obj:`bids.layout.models.Config` objects or names of
        configurations to use in matching. Each element must be a
        :obj:`bids.layout.models.Config` object, or a valid
        :obj:`bids.layout.models.Config` name (e.g., 'bids' or 'derivatives').
        If None, all available configs are used. Default is None.
    include_unmatched : bool, optional
        If True, unmatched entities are included in the returned dict,
        with values set to None.
        If False (default), unmatched entities are ignored.

    Returns
    -------
    dict
        Keys are Entity names and values are the values from the filename.
    """
    # Load Configs if needed
    if entities is None:

        if config is None:
            config = ['bids', 'derivatives']

        config = [Config.load(c) if not isinstance(c, Config) else c
                  for c in listify(config)]

        # Consolidate entities from all Configs into a single dict
        entities = {}
        for c in config:
            entities.update(c.entities)
        entities = entities.values()

    # Extract matches
    bf = make_bidsfile(filename)
    ent_vals = {}
    for ent in entities:
        match = ent.match_file(bf)
        if match is not None or include_unmatched:
            ent_vals[ent.name] = match

    return ent_vals


def add_config_paths(**kwargs):
    """Add to the pool of available configuration files for BIDSLayout.

    Parameters
    ----------
    kwargs : dict
        Dictionary specifying where to find additional config files.
        Keys are names, values are paths to the corresponding .json file.

    Examples
    --------
    > add_config_paths(my_config='/path/to/config')
    > layout = BIDSLayout('/path/to/bids', config=['bids', 'my_config'])
    """
    for k, path in kwargs.items():
        if not os.path.exists(path):
            raise ValueError(
                'Configuration file "{}" does not exist'.format(k))
        if k in cf.get_option('config_paths'):
            raise ValueError('Configuration {!r} already exists'.format(k))

    kwargs.update(**cf.get_option('config_paths'))
    cf.set_option('config_paths', kwargs)


class BIDSLayout(object):
    """Layout class representing an entire BIDS dataset.

    Parameters
    ----------
    root : str
        The root directory of the BIDS dataset.
    validate : bool, optional
        If True, all files are checked for BIDS compliance
        when first indexed, and non-compliant files are ignored. This
        provides a convenient way to restrict file indexing to only those
        files defined in the "core" BIDS spec, as setting validate=True
        will lead files in supplementary folders like derivatives/, code/,
        etc. to be ignored.
    absolute_paths : bool, optional
        If True, queries always return absolute paths.
        If False, queries return relative paths (for files and
        directories).
    derivatives : bool or str or list, optional
        Specifies whether and/or which
        derivatives to to index. If True, all pipelines found in the
        derivatives/ subdirectory will be indexed. If a str or list, gives
        the paths to one or more derivatives directories to index. If False
        or None, the derivatives/ directory is ignored during indexing, and
        derivatives will have to be added manually via add_derivatives().
        Note: derivatives datasets MUST contain a dataset_description.json
        file in order to be indexed.
    config : str or list or None, optional
        Optional name(s) of configuration file(s) to use.
        By default (None), uses 'bids'.
    sources : :obj:`bids.layout.BIDSLayout` or list or None, optional
        Optional BIDSLayout(s) from which the current BIDSLayout is derived.
    ignore : str or SRE_Pattern or list
        Path(s) to exclude from indexing. Each
        path is either a string or a SRE_Pattern object (i.e., compiled
        regular expression). If a string is passed, it must be either an
        absolute path, or be relative to the BIDS project root. If an
        SRE_Pattern is passed, the contained regular expression will be
        matched against the full (absolute) path of all files and
        directories. By default, indexing ignores all files in 'code/',
        'stimuli/', 'sourcedata/', 'models/', and any hidden files/dirs
        beginning with '.' at root level.
    force_index : str or SRE_Pattern or list
        Path(s) to forcibly index in the
        BIDSLayout, even if they would otherwise fail validation. See the
        documentation for the ignore argument for input format details.
        Note that paths in force_index takes precedence over those in
        ignore (i.e., if a file matches both ignore and force_index, it
        *will* be indexed).
        Note: NEVER include 'derivatives' here; use the derivatives argument
        (or :obj:`bids.layout.BIDSLayout.add_derivatives`) for that.
    config_filename : str
        Optional name of filename within directories
        that contains configuration information.
    regex_search : bool
        Whether to require exact matching (True) or regex
        search (False, default) when comparing the query string to each
        entity in .get() calls. This sets a default for the instance, but
        can be overridden in individual .get() requests.
    database_path : str
        Optional path to directory containing SQLite database file index
        for this BIDS dataset. If a value is passed and the folder
        already exists, indexing is skipped. By default (i.e., if None),
        an in-memory SQLite database is used, and the index will not
        persist unless .save() is explicitly called.
    reset_database : bool
        If True, any existing directory specified in the
        database_path argument is deleted, and the BIDS dataset provided
        in the root argument is reindexed. If False, indexing will be
        skipped and the existing database file will be used. Ignored if
        database_path is not provided.
    index_metadata : bool
        If True, all metadata files are indexed at
        initialization. If False, metadata will not be available (but
        indexing will be faster).
    """

    _default_ignore = ("code", "stimuli", "sourcedata", "models",
                       re.compile(r'^\.'))

    def __init__(self, root, validate=False, absolute_paths=True,
                 derivatives=False, config=None, sources=None, ignore=None,
                 force_index=None, config_filename='layout_config.json',
                 regex_search=False, database_path=None, database_file=None,
                 reset_database=False, index_metadata=True):
        """Initialize BIDSLayout."""
        self.root = root
        self.validate = validate
        self.absolute_paths = absolute_paths
        self.derivatives = {}
        self.sources = sources
        self.regex_search = regex_search
        self.config_filename = config_filename
        # Store original init arguments as dictionary
        self._init_args = self._sanitize_init_args(
            root=root, validate=validate, absolute_paths=absolute_paths,
            derivatives=derivatives, ignore=ignore, force_index=force_index,
            index_metadata=index_metadata, config=config)

        if database_path is None and database_file is not None:
            database_path = database_file
            warnings.warn(
                'In pybids 0.10 database_file argument was deprecated in favor'
                ' of database_path, and will be removed in 0.12. '
                'For now, treating database_file as a directory.',
                DeprecationWarning)
        if database_path:
            database_path = str(Path(database_path).absolute())

        self.session = None

        index_dataset = self._init_db(database_path, reset_database)

        # Do basic BIDS validation on root directory
        self._validate_root()

        if ignore is None:
            ignore = self._default_ignore

        # Instantiate after root validation to ensure os.path.join works
        self.ignore = [os.path.abspath(os.path.join(self.root, patt))
                       if isinstance(patt, str) else patt
                       for patt in listify(ignore or [])]
        self.force_index = [os.path.abspath(os.path.join(self.root, patt))
                            if isinstance(patt, str) else patt
                            for patt in listify(force_index or [])]

        # Initialize the BIDS validator and examine ignore/force_index args
        self._validate_force_index()

        if index_dataset:
            # Create Config objects
            if config is None:
                config = 'bids'
            config = [Config.load(c, session=self.session)
                      for c in listify(config)]
            self.config = {c.name: c for c in config}
            # Missing persistence of configs to the database
            for config_obj in self.config.values():
                self.session.add(config_obj)
                self.session.commit()

            # Index files and (optionally) metadata
            indexer = BIDSLayoutIndexer(self)
            indexer.index_files()
            if index_metadata:
                indexer.index_metadata()
        else:
            # Load Configs from DB
            self.config = {c.name: c for c in self.session.query(Config).all()}

        # Add derivatives if any are found
        if derivatives:
            if derivatives is True:
                derivatives = os.path.join(root, 'derivatives')
            self.add_derivatives(
                derivatives, parent_database_path=database_path,
                validate=validate, absolute_paths=absolute_paths,
                derivatives=None, sources=self, ignore=ignore,  config=None,
                force_index=force_index, config_filename=config_filename,
                regex_search=regex_search, index_metadata=index_metadata,
                reset_database=index_dataset or reset_database
                )

    def __getattr__(self, key):
        """Dynamically inspect missing methods for get_<entity>() calls
        and return a partial function of get() if a match is found."""
        if key.startswith('get_'):
            ent_name = key.replace('get_', '')
            entities = self.get_entities()
            # Use inflect to check both singular and plural forms
            if ent_name not in entities:
                sing = inflect.engine().singular_noun(ent_name)
                if sing in entities:
                    ent_name = sing
                else:
                    raise AttributeError(
                        "'get_{}' can't be called because '{}' isn't a "
                        "recognized entity name.".format(ent_name, ent_name))
            return partial(self.get, return_type='id', target=ent_name)
        # Spit out default message if we get this far
        raise AttributeError("%s object has no attribute named %r" %
                             (self.__class__.__name__, key))

    def __repr__(self):
        """Provide a tidy summary of key properties."""
        n_subjects = len(
            [s.value
             for s in self.session.query(Tag).filter_by(
                 entity_name='subject').group_by(Tag._value)]
            )

        n_sessions = len(
            set(
                (t.value, t.file.entities.get('subject'))
                 for t in
                 self.session.query(Tag).filter_by(entity_name='session')
                 if t.file.entities.get('subject')
                 )
            )

        n_runs = len(
            set(
                (t.value, t.file.entities.get('subject'))
                 for t in
                 self.session.query(Tag).filter_by(entity_name='run')
                 if t.file.entities.get('subject')
                 )
            )

        root = self.root[-30:]
        s = ("BIDS Layout: ...{} | Subjects: {} | Sessions: {} | "
             "Runs: {}".format(root, n_subjects, n_sessions, n_runs))
        return s

    def _set_session(self, database_file):
        if database_file is not None:
            # https://docs.sqlalchemy.org/en/13/dialects/sqlite.html
            # When a file-based database is specified, the dialect will use
            # NullPool as the source of connections. This pool closes and
            # discards connections which are returned to the pool immediately.
            # SQLite file-based connections have extremely low overhead, so
            # pooling is not necessary. The scheme also prevents a connection
            # from being used again in a different thread and works best
            # with SQLite's coarse-grained file locking.
            from sqlalchemy.pool import NullPool
            engine = sa.create_engine(
                'sqlite:///{dbfilepath}'.format(dbfilepath=database_file),
                connect_args={'check_same_thread': False},
                poolclass=NullPool)
        else:
            # https://docs.sqlalchemy.org/en/13/dialects/sqlite.html
            # Using a Memory Database in Multiple Threads
            # To use a :memory: database in a multithreaded scenario, the same
            # connection object must be shared among
            # threads, since the database exists only within the scope of that
            # connection. The StaticPool implementation will maintain a single
            # connection globally, and the check_same_thread flag can be passed
            # to Pysqlite as False:
            from sqlalchemy.pool import StaticPool
            engine = sa.create_engine(
                'sqlite://',  # In memory database
                connect_args={'check_same_thread': False},
                poolclass=StaticPool)
            # Note that using a :memory: database in multiple threads requires
            # a recent version of SQLite.

        def regexp(expr, item):
            """Regex function for SQLite's REGEXP."""
            reg = re.compile(expr, re.I)
            return reg.search(item) is not None

        conn = engine.connect()

        # Do not remove this decorator!!! An in-line create_function call will
        # work when using an in-memory SQLite DB, but fails when using a file.
        # For more details, see https://stackoverflow.com/questions/12461814/
        @sa.event.listens_for(engine, "begin")
        def do_begin(conn):
            conn.connection.create_function('regexp', 2, regexp)

        self.session = sa.orm.sessionmaker(bind=engine)()

    @staticmethod
    def _make_db_paths(database_path):
        if database_path is not None:
            database_path = Path(database_path)
            database_file = database_path / 'layout_index.sqlite'
            database_sidecar = database_path / 'layout_args.json'
            database_path.mkdir(exist_ok=True, parents=True)
        else:
            database_file = None
            database_sidecar = None
        return database_file, database_sidecar

    @staticmethod
    def _sanitize_init_args(**kwargs):
        """ Prepare initalization arguments for serialization """
        # Make ignore and force_index serializable
        for k in ['ignore', 'force_index']:
            kwargs[k] = [
                str(a) for a in kwargs[k] if a is not None] \
              if kwargs[k] is not None else None

        kwargs['root'] = str(Path(kwargs['root']).absolute())

        # Get abspaths
        if isinstance(kwargs['derivatives'], list):
            kwargs['derivatives'] = [
                str(Path(der).absolute())
                for der in listify(kwargs['derivatives'])
                ]

        return kwargs

    def _init_db(self, database_path=None, reset_database=False):
        database_file, database_sidecar = self._make_db_paths(database_path)
        # Reset database if needed and return whether or not it was reset.
        # Determining if the database needs resetting must be done prior
        # to setting the session (which creates the empty database file)
        reset_database = (
            reset_database or  # Manual Request
            not database_file or  # In memory transient db
            not database_file.exists()  # New file based db created
        )

        self._set_session(database_file)

        if not reset_database:
            saved_args = json.loads(database_sidecar.read_text())
            for k, v in saved_args.items():
                if self._init_args[k] != v:
                    raise ValueError(
                        "Initialization argument ('{}') does not match "
                        "for database_path: {}.\n"
                        "Saved value: {}.\n"
                        "Current value: {}.".format(
                            k, database_path, v, self._init_args[k])
                        )
        else:
            engine = self.session.get_bind()
            Base.metadata.drop_all(engine)
            Base.metadata.create_all(engine)
            if database_sidecar:
                database_sidecar.write_text(json.dumps(self._init_args))

            return True

        return False

    def _validate_root(self):
        # Validate root argument and make sure it contains mandatory info

        try:
            self.root = str(self.root)
        except:
            raise TypeError("root argument must be a string (or a type that "
                            "supports casting to string, such as "
                            "pathlib.Path) specifying the directory "
                            "containing the BIDS dataset.")

        self.root = os.path.abspath(self.root)

        if not os.path.exists(self.root):
            raise ValueError("BIDS root does not exist: %s" % self.root)

        target = os.path.join(self.root, 'dataset_description.json')
        if not os.path.exists(target):
            if self.validate:
                raise ValueError(
                    "'dataset_description.json' is missing from project root."
                    " Every valid BIDS dataset must have this file.")
            else:
                self.description = None
        else:
            with open(target, 'r', encoding='utf-8') as desc_fd:
                self.description = json.load(desc_fd)
            if self.validate:
                for k in ['Name', 'BIDSVersion']:
                    if k not in self.description:
                        raise ValueError("Mandatory '%s' field missing from "
                                         "dataset_description.json." % k)

    def _validate_force_index(self):
        # Derivatives get special handling; they shouldn't be indexed normally
        if self.force_index is not None:
            for entry in self.force_index:
                condi = (isinstance(entry, str) and
                         os.path.normpath(entry).startswith('derivatives'))
                if condi:
                    msg = ("Do not pass 'derivatives' in the force_index "
                           "list. To index derivatives, either set "
                           "derivatives=True, or use add_derivatives().")
                    raise ValueError(msg)

    def _in_scope(self, scope):
        """Determine whether current BIDSLayout is in the passed scope.

        Parameters
        ----------
        scope : str or list
            The intended scope(s). Each value must be one of 'all', 'raw',
            'derivatives', or a pipeline name.
        """
        scope = listify(scope)

        if 'all' in scope:
            return True

        # We assume something is a BIDS-derivatives dataset if it either has a
        # defined pipeline name, or is applying the 'derivatives' rules.
        pl_name = self.description.get("PipelineDescription", {}).get("Name")
        is_deriv = bool('derivatives' in self.config)

        return ((not is_deriv and 'raw' in scope) or
                (is_deriv and ('derivatives' in scope or pl_name in scope)))

    def _get_layouts_in_scope(self, scope):
        """Return all layouts in the passed scope."""

        if scope == 'self':
            return [self]

        def collect_layouts(layout):
            """Recursively build a list of layouts."""
            children = list(layout.derivatives.values())
            layouts = [collect_layouts(d) for d in children]
            return [layout] + list(chain(*layouts))

        layouts = [l for l in collect_layouts(self) if l._in_scope(scope)]
        return list(set(layouts))

    def _sanitize_query_dtypes(self, entities):
        """Automatically convert entity query values to correct dtypes."""
        entities = entities.copy()
        names = list(entities.keys())
        ents = {e.name: e for e in
                self.session.query(Entity)
                    .filter(Entity.name.in_(names)).all()}
        # Fail silently because the DB may still know how to reconcile
        # type differences.
        for name, val in entities.items():
            if isinstance(val, enum.Enum):
                continue
            try:
                if isinstance(val, (list, tuple)):
                    entities[name] = [ents[name]._astype(v) for v in val]
                else:
                    entities[name] = ents[name]._astype(val)
            except:
                pass
        return entities

    @property
    def entities(self):
        """Get the entities."""
        return self.get_entities()

    @property
    def files(self):
        """Get the files."""
        return self.get_files()

    @classmethod
    def load(cls, database_path):
        """ Load index from database path. Initalization parameters are set to
        those found in database_path JSON sidecar.

        Parameters
        ----------
        database_path : str
            The path to the desired database folder. By default,
            uses .db_cache. If a relative path is passed, it is assumed to
            be relative to the BIDSLayout root directory.
        """
        database_file, database_sidecar = cls._make_db_paths(database_path)
        init_args = json.loads(database_sidecar.read_text())

        return cls(database_path=database_path, **init_args)

    def save(self, database_path, replace_connection=True):
        """Save the current index as a SQLite3 DB at the specified location.

        Note: This is only necessary if a database_path was not specified
        at initialization, and the user now wants to save the index.
        If a database_path was specified originally, there is no need to
        re-save using this method.

        Parameters
        ----------
        database_path : str
            The path to the desired database folder. By default,
            uses .db_cache. If a relative path is passed, it is assumed to
            be relative to the BIDSLayout root directory.
        replace_connection : bool, optional
            If True, the newly created database will
            be used for all subsequent connections. This means that any
            changes to the index made after the .save() call will be
            reflected in the database file. If False, the previous database
            will continue to be used, and any subsequent changes will not
            be reflected in the new file unless save() is explicitly called
            again.
        """
        database_file, database_sidecar = self._make_db_paths(database_path)
        new_db = sqlite3.connect(str(database_file))
        old_db = self.session.get_bind().connect().connection

        with new_db:
            for line in old_db.iterdump():
                if line not in ('BEGIN;', 'COMMIT;'):
                    new_db.execute(line)
            new_db.commit()

        if replace_connection:
            self._set_session(str(database_file))

        # Dump instance arguments to JSON
        database_sidecar.write_text(json.dumps(self._init_args))

        # Recursively save children
        for pipeline_name, der in self.derivatives.items():
            der.save(os.path.join(
                database_path, pipeline_name))

    def get_entities(self, scope='all', metadata=None):
        """Get entities for all layouts in the specified scope.

        Parameters
        ----------
        scope : str
            The scope of the search space. Indicates which
            BIDSLayouts' entities to extract.
            See :obj:`bids.layout.BIDSLayout.get` docstring for valid values.

        metadata : bool or None
            By default (None), all available entities
            are returned. If True, only entities found in metadata files
            (and not defined for filenames) are returned. If False, only
            entities defined for filenames (and not those found in JSON
            sidecars) are returned.

        Returns
        -------
        dict
            Dictionary where keys are entity names and
            values are Entity instances.
        """
        # TODO: memoize results
        layouts = self._get_layouts_in_scope(scope)
        entities = {}
        for l in layouts:
            query = l.session.query(Entity)
            if metadata is not None:
                query = query.filter_by(is_metadata=metadata)
            results = query.all()
            entities.update({e.name: e for e in results})
        return entities

    def get_files(self, scope='all'):
        """Get BIDSFiles for all layouts in the specified scope.

        Parameters
        ----------
        scope : str
            The scope of the search space. Indicates which
            BIDSLayouts' entities to extract.
            See :obj:`bids.layout.BIDSLayout.get` docstring for valid values.


        Returns:
            A dict, where keys are file paths and values
            are :obj:`bids.layout.BIDSFile` instances.

        """
        # TODO: memoize results
        layouts = self._get_layouts_in_scope(scope)
        files = {}
        for l in layouts:
            results = l.session.query(BIDSFile).all()
            files.update({f.path: f for f in results})
        return files

    def clone(self):
        """Return a deep copy of the current BIDSLayout."""
        return copy.deepcopy(self)

    def parse_file_entities(self, filename, scope='all', entities=None,
                            config=None, include_unmatched=False):
        """Parse the passed filename for entity/value pairs.

        Parameters
        ----------
        filename : str
            The filename to parse for entity values
        scope : str or list, optional
            The scope of the search space. Indicates which BIDSLayouts'
            entities to extract. See :obj:`bids.layout.BIDSLayout.get`
            docstring for valid values. By default, extracts all entities.
        entities : list or None, optional
            An optional list of Entity instances to use in
            extraction. If passed, the scope and config arguments are
            ignored, and only the Entities in this list are used.
        config : str or :obj:`bids.layout.models.Config` or list or None, optional
            One or more :obj:`bids.layout.models.Config` objects, or paths
            to JSON config files on disk, containing the Entity definitions
            to use in extraction. If passed, scope is ignored.
        include_unmatched : bool, optional
            If True, unmatched entities are included
            in the returned dict, with values set to None. If False
            (default), unmatched entities are ignored.

        Returns
        -------
        dict
            Dictionary where keys are Entity names and values are the
            values extracted from the filename.
        """
        # If either entities or config is specified, just pass through
        if entities is None and config is None:
            layouts = self._get_layouts_in_scope(scope)
            config = chain(*[list(l.config.values()) for l in layouts])
            config = list(set(config))

        return parse_file_entities(filename, entities, config,
                                   include_unmatched)

    def add_derivatives(self, path, parent_database_path=None, **kwargs):
        """Add BIDS-Derivatives datasets to tracking.

        Parameters
        ----------
        path : str or list
            One or more paths to BIDS-Derivatives datasets.
            Each path can point to either a derivatives/ directory
            containing one more more pipeline directories, or to a single
            pipeline directory (e.g., derivatives/fmriprep).
        parent_database_path : str
            If not None, use the pipeline name from the dataset_description.json
            file as the database folder name to nest within the parent database
            folder name to write out derivative index to.
        kwargs : dict
            Optional keyword arguments to pass on to
            BIDSLayout() when initializing each of the derivative datasets.

        Notes
        -----
        Every derivatives directory intended for indexing MUST contain a
        valid dataset_description.json file. See the BIDS-Derivatives
        specification for details.
        """
        paths = listify(path)
        deriv_dirs = []

        # Collect all paths that contain a dataset_description.json
        def check_for_description(bids_dir):
            dd = os.path.join(bids_dir, 'dataset_description.json')
            return os.path.exists(dd)

        for p in paths:
            p = os.path.abspath(p)
            if os.path.exists(p):
                if check_for_description(p):
                    deriv_dirs.append(p)
                else:
                    subdirs = [d for d in os.listdir(p)
                               if os.path.isdir(os.path.join(p, d))]
                    for sd in subdirs:
                        sd = os.path.join(p, sd)
                        if check_for_description(sd):
                            deriv_dirs.append(sd)

        if not deriv_dirs:
            warnings.warn("Derivative indexing was requested, but no valid "
                          "datasets were found in the specified locations "
                          "({}). Note that all BIDS-Derivatives datasets must"
                          " meet all the requirements for BIDS-Raw datasets "
                          "(a common problem is to fail to include a "
                          "dataset_description.json file in derivatives "
                          "datasets).".format(paths))

        for deriv in deriv_dirs:
            dd = os.path.join(deriv, 'dataset_description.json')
            with open(dd, 'r', encoding='utf-8') as ddfd:
                description = json.load(ddfd)
            pipeline_name = description.get(
                'PipelineDescription', {}).get('Name')
            if pipeline_name is None:
                raise ValueError("Every valid BIDS-derivatives dataset must "
                                 "have a PipelineDescription.Name field set "
                                 "inside dataset_description.json.")
            if pipeline_name in self.derivatives:
                raise ValueError("Pipeline name '%s' has already been added "
                                 "to this BIDSLayout. Every added pipeline "
                                 "must have a unique name!")
            # Default config and sources values
            kwargs['config'] = kwargs.get('config') or ['bids', 'derivatives']
            kwargs['sources'] = kwargs.get('sources') or self
            if parent_database_path:
                child_database_path = os.path.join(
                    parent_database_path, pipeline_name)
                kwargs['database_path'] = child_database_path

            self.derivatives[pipeline_name] = BIDSLayout(deriv, **kwargs)

    def to_df(self, metadata=False, **filters):
        """Return information for BIDSFiles tracked in Layout as pd.DataFrame.

        Parameters
        ----------
        metadata : bool, optional
            If True, includes columns for all metadata fields.
            If False, only filename-based entities are included as columns.
        filters : dict, optional
            Optional keyword arguments passed on to get(). This allows
            one to easily select only a subset of files for export.

        Returns
        -------
        :obj:`pandas.DataFrame`
            A pandas DataFrame, where each row is a file, and each column is
            a tracked entity. NaNs are injected whenever a file has no
            value for a given attribute.

        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError('Missing dependency: "pandas"')

        # TODO: efficiency could probably be improved further by joining the
        # BIDSFile and Tag tables and running a single query. But this would
        # require refactoring the below to use _build_file_query, which will
        # in turn likely require generalizing the latter.
        files = self.get(**filters)
        file_paths = [f.path for f in files]
        query = self.session.query(Tag).filter(Tag.file_path.in_(file_paths))

        if not metadata:
            query = query.join(Entity).filter(Entity.is_metadata == False)

        tags = query.all()

        tags = [[t.file_path, t.entity_name, t.value] for t in tags]
        data = pd.DataFrame(tags, columns=['path', 'entity', 'value'])
        data = data.pivot('path', 'entity', 'value')

        # Add in orphaned files with no Tags. Maybe make this an argument?
        orphans = list(set(file_paths) - set(data.index))
        for o in orphans:
            data.loc[o] = pd.Series()

        return data.reset_index()

    def get(self, return_type='object', target=None, scope='all',
            regex_search=False, absolute_paths=None, drop_invalid_filters=True,
            **filters):
        """Retrieve files and/or metadata from the current Layout.

        Parameters
        ----------
        return_type : str, optional
            Type of result to return. Valid values:
            'object' (default): return a list of matching BIDSFile objects.
            'file' or 'filename': return a list of matching filenames.
            'dir': return a list of directories.
            'id': return a list of unique IDs. Must be used together
                  with a valid target.
        target : str, optional
            Optional name of the target entity to get results for
            (only used if return_type is 'dir' or 'id').
        scope : str or list, optional
            Scope of the search space. If passed, only
            nodes/directories that match the specified scope will be
            searched. Possible values include:
            'all' (default): search all available directories.
            'derivatives': search all derivatives directories.
            'raw': search only BIDS-Raw directories.
            'self': search only the directly called BIDSLayout.
            <PipelineName>: the name of a BIDS-Derivatives pipeline.
        regex_search : bool or None, optional
            Whether to require exact matching
            (False) or regex search (True) when comparing the query string
            to each entity.
        absolute_paths : bool, optional
            Optionally override the instance-wide option
            to report either absolute or relative (to the top of the
            dataset) paths. If None, will fall back on the value specified
            at BIDSLayout initialization.
        filters : dict
            Any optional key/values to filter the entities on.
            Keys are entity names, values are regexes to filter on. For
            example, passing filters={'subject': 'sub-[12]'} would return
            only files that match the first two subjects. In addition to
            ordinary data types, the following enums are defined (in the
            Query class):
                * Query.NONE: The named entity must not be defined.
                * Query.ANY: the named entity must be defined, but can have any
                    value.

        Returns
        -------
        list of :obj:`bids.layout.BIDSFile` or str
            A list of BIDSFiles (default) or strings (see return_type).

        Notes
        -----
        * In pybids 0.9.0, the 'extensions' argument has been removed in
          favor of the 'extension' entity.
        """
        if 'extensions' in filters:
            filters['extension'] = filters.pop('extensions')
            warnings.warn("In pybids 0.9.0, the 'extensions' filter was "
                          "deprecated in favor of 'extension'. The former will"
                          " stop working in 0.11.0.", DeprecationWarning)

        layouts = self._get_layouts_in_scope(scope)

        entities = self.get_entities()

        # For consistency with past versions where "extensions" was a
        # hard-coded argument, allow leading periods
        if 'extension' in filters:
            exts = listify(filters['extension'])
            filters['extension'] = [x.lstrip('.') for x in exts]

        if drop_invalid_filters:
            invalid_filters = set(filters.keys()) - set(entities.keys())
            if invalid_filters:
                for inv_filt in invalid_filters:
                    filters.pop(inv_filt)

        # Provide some suggestions if target is specified and invalid.
        if target is not None and target not in entities:
            import difflib
            potential = list(entities.keys())
            suggestions = difflib.get_close_matches(target, potential)
            if suggestions:
                message = "Did you mean one of: {}?".format(suggestions)
            else:
                message = "Valid targets are: {}".format(potential)
            raise ValueError(("Unknown target '{}'. " + message)
                             .format(target))

        results = []
        for l in layouts:
            query = l._build_file_query(filters=filters,
                                        regex_search=regex_search)
            # Eager load associations, because mixing queries from different
            # DB sessions causes objects to detach
            query = query.options(joinedload(BIDSFile.tags)
                                  .joinedload(Tag.entity))
            results.extend(query.all())

        # Convert to relative paths if needed
        if absolute_paths is None:  # can be overloaded as option to .get
            absolute_paths = self.absolute_paths

        if not absolute_paths:
            for i, fi in enumerate(results):
                fi = copy.copy(fi)
                fi.path = os.path.relpath(fi.path, self.root)
                results[i] = fi

        if return_type.startswith('file'):
            results = natural_sort([f.path for f in results])

        elif return_type in ['id', 'dir']:
            if target is None:
                raise ValueError('If return_type is "id" or "dir", a valid '
                                 'target entity must also be specified.')

            results = [x for x in results if target in x.entities]

            if return_type == 'id':
                results = list(set([x.entities[target] for x in results]))
                results = natural_sort(results)

            elif return_type == 'dir':
                template = entities[target].directory
                if template is None:
                    raise ValueError('Return type set to directory, but no '
                                     'directory template is defined for the '
                                     'target entity (\"%s\").' % target)
                # Construct regex search pattern from target directory template
                template = self.root + template
                to_rep = re.findall(r'{(.*?)\}', template)
                for ent in to_rep:
                    patt = entities[ent].pattern
                    template = template.replace('{%s}' % ent, patt)
                template += r'[^\%s]*$' % os.path.sep
                matches = [
                    f.dirname if absolute_paths else os.path.relpath(f.dirname, self.root)  # noqa: E501
                    for f in results
                    if re.search(template, f.dirname)
                ]

                results = natural_sort(list(set(matches)))

            else:
                raise ValueError("Invalid return_type specified (must be one "
                                 "of 'tuple', 'filename', 'id', or 'dir'.")
        else:
            results = natural_sort(results, 'path')

        return results

    def get_file(self, filename, scope='all'):
        """Return the BIDSFile object with the specified path.

        Parameters
        ----------
        filename : str
            The path of the file to retrieve. Must be either an absolute path,
            or relative to the root of this BIDSLayout.
        scope : str or list, optional
            Scope of the search space. If passed, only BIDSLayouts that match
            the specified scope will be searched. See :obj:`BIDSLayout.get`
            docstring for valid values. Default is 'all'.

        Returns
        -------
        :obj:`bids.layout.BIDSFile` or None
            File found, or None if no match was found.
        """
        filename = os.path.abspath(os.path.join(self.root, filename))
        for layout in self._get_layouts_in_scope(scope):
            result = layout.session.query(
                BIDSFile).filter_by(path=filename).first()  # noqa: E501
            if result:
                return result
        return None

    def _build_file_query(self, **kwargs):

        query = self.session.query(BIDSFile).filter_by(is_dir=False)

        filters = kwargs.get('filters')

        # Entity filtering
        if filters:
            query = query.join(BIDSFile.tags)
            regex = kwargs.get('regex_search', False)

            filters = self._sanitize_query_dtypes(filters)

            for name, val in filters.items():
                if isinstance(val, (list, tuple)) and len(val) == 1:
                    val = val[0]
                if val is None or isinstance(val, enum.Enum):
                    name_clause = query.filter(BIDSFile.tags.any(entity_name=name))
                    if val == Query.ANY:
                        query = name_clause
                    else:
                        query = query.except_(name_clause)
                    continue

                if regex:
                    if isinstance(val, (list, tuple)):
                        val_clause = sa.or_(*[Tag._value.op('REGEXP')(str(v))
                                              for v in val])
                    else:
                        val_clause = Tag._value.op('REGEXP')(str(val))
                else:
                    if isinstance(val, (list, tuple)):
                        val_clause = Tag._value.in_(val)
                    else:
                        val_clause = Tag._value == val

                subq = sa.and_(Tag.entity_name == name, val_clause)
                query = query.filter(BIDSFile.tags.any(subq))

        return query

    def get_collections(self, level, types=None, variables=None, merge=False,
                        sampling_rate=None, skip_empty=False, **kwargs):
        """Return one or more variable Collections in the BIDS project.

        Parameters
        ----------
        level : {'run', 'session', 'subject', 'dataset'}
            The level of analysis to return variables for.
            Must be one of 'run', 'session','subject', or 'dataset'.
        types : str or list
            Types of variables to retrieve. All valid values reflect the
            filename stipulated in the BIDS spec for each kind of variable.
            Valid values include: 'events', 'physio', 'stim', 'scans',
            'participants', 'sessions', and 'regressors'. Default is None.
        variables : list
            Optional list of variables names to return. If None, all available
            variables are returned.
        merge : bool
            If True, variables are merged across all observations of the
            current level. E.g., if level='subject', variables from all
            subjects will be merged into a single collection. If False, each
            observation is handled separately, and the result is returned
            as a list.
        sampling_rate : int or str
            If level='run', the sampling rate to pass on to the returned
            :obj:`bids.variables.kollekshuns.BIDSRunVariableCollection`.
        skip_empty : bool
            Whether or not to skip empty Variables (i.e., where there are no
            rows/records in a file after applying any filtering operations
            like dropping NaNs).
        kwargs
            Optional additional arguments to pass onto
            :obj:`bids.variables.io.load_variables`.

        Returns
        -------
        list of :obj:`bids.variables.kollekshuns.BIDSVariableCollection`
            or :obj:`bids.variables.kollekshuns.BIDSVariableCollection`
            A list if merge=False;
            a single :obj:`bids.variables.kollekshuns.BIDSVariableCollection`
            if merge=True.

        """
        from bids.variables import load_variables
        index = load_variables(self, types=types, levels=level,
                               skip_empty=skip_empty, **kwargs)
        return index.get_collections(level, variables, merge,
                                     sampling_rate=sampling_rate)

    def get_metadata(self, path, include_entities=False, scope='all'):
        """Return metadata found in JSON sidecars for the specified file.

        Parameters
        ----------
        path : str
            Path to the file to get metadata for.
        include_entities : bool, optional
            If True, all available entities extracted
            from the filename (rather than JSON sidecars) are included in
            the returned metadata dictionary.
        scope : str or list, optional
            The scope of the search space. Each element must
            be one of 'all', 'raw', 'self', 'derivatives', or a
            BIDS-Derivatives pipeline name. Defaults to searching all
            available datasets.

        Returns
        -------
        dict
            A dictionary of key/value pairs extracted from all of the
            target file's associated JSON sidecars.

        Notes
        -----
        A dictionary containing metadata extracted from all matching .json
        files is returned. In cases where the same key is found in multiple
        files, the values in files closer to the input filename will take
        precedence, per the inheritance rules in the BIDS specification.

        """
        for layout in self._get_layouts_in_scope(scope):

            query = (layout.session.query(Tag)
                     .join(BIDSFile)
                     .filter(BIDSFile.path == path))

            if not include_entities:
                query = query.join(Entity).filter(Entity.is_metadata == True)

            results = query.all()
            if results:
                return {t.entity_name: t.value for t in results}

        return {}

    def get_dataset_description(self, scope='self', all_=False):
        """Return contents of dataset_description.json.

        Parameters
        ----------
        scope : str
            The scope of the search space. Only descriptions of
            BIDSLayouts that match the specified scope will be returned.
            See :obj:`bids.layout.BIDSLayout.get` docstring for valid values.
            Defaults to 'self' --i.e., returns the dataset_description.json
            file for only the directly-called BIDSLayout.
        all_ : bool
            If True, returns a list containing descriptions for
            all matching layouts. If False (default), returns for only the
            first matching layout.

        Returns
        -------
        dict or list of dict
            a dictionary or list of dictionaries (depending on all_).
        """
        layouts = self._get_layouts_in_scope(scope)
        if not all_:
            return layouts[0].get_file('dataset_description.json').get_dict()
        return [l.get_file('dataset_description.json').get_dict()
                for l in layouts]

    def get_nearest(self, path, return_type='filename', strict=True,
                    all_=False, ignore_strict_entities='extension',
                    full_search=False, **filters):
        """Walk up file tree from specified path and return nearest matching file(s).

        Parameters
        ----------
        path (str): The file to search from.
        return_type (str): What to return; must be one of 'filename'
            (default) or 'tuple'.
        strict (bool): When True, all entities present in both the input
            path and the target file(s) must match perfectly. When False,
            files will be ordered by the number of matching entities, and
            partial matches will be allowed.
        all_ (bool): When True, returns all matching files. When False
            (default), only returns the first match.
        ignore_strict_entities (str, list): Optional entity/entities to
            exclude from strict matching when strict is True. This allows
            one to search, e.g., for files of a different type while
            matching all other entities perfectly by passing
            ignore_strict_entities=['type']. Ignores extension by default.
        full_search (bool): If True, searches all indexed files, even if
            they don't share a common root with the provided path. If
            False, only files that share a common root will be scanned.
        filters : dict
            Optional keywords to pass on to :obj:`bids.layout.BIDSLayout.get`.
        """
        path = os.path.abspath(path)

        # Make sure we have a valid suffix
        if not filters.get('suffix'):
            f = self.get_file(path)
            if 'suffix' not in f.entities:
                raise ValueError(
                    "File '%s' does not have a valid suffix, most "
                    "likely because it is not a valid BIDS file." % path
                )
            filters['suffix'] = f.entities['suffix']

        # Collect matches for all entities
        entities = {}
        for ent in self.get_entities(metadata=False).values():
            m = ent.regex.search(path)
            if m:
                entities[ent.name] = ent._astype(m.group(1))

        # Remove any entities we want to ignore when strict matching is on
        if strict and ignore_strict_entities is not None:
            for k in listify(ignore_strict_entities):
                entities.pop(k, None)

        # Get candidate files
        results = self.get(**filters)

        # Make a dictionary of directories --> contained files
        folders = defaultdict(list)
        for f in results:
            folders[f.dirname].append(f)

        # Build list of candidate directories to check
        search_paths = []
        while True:
            if path in folders and folders[path]:
                search_paths.append(path)
            parent = os.path.dirname(path)
            if parent == path:
                break
            path = parent

        if full_search:
            unchecked = set(folders.keys()) - set(search_paths)
            search_paths.extend(path for path in unchecked if folders[path])

        def count_matches(f):
            # Count the number of entities shared with the passed file
            f_ents = f.entities
            keys = set(entities.keys()) & set(f_ents.keys())
            shared = len(keys)
            return [shared, sum([entities[k] == f_ents[k] for k in keys])]

        matches = []

        for path in search_paths:
            # Sort by number of matching entities. Also store number of
            # common entities, for filtering when strict=True.
            num_ents = [[f] + count_matches(f) for f in folders[path]]
            # Filter out imperfect matches (i.e., where number of common
            # entities does not equal number of matching entities).
            if strict:
                num_ents = [f for f in num_ents if f[1] == f[2]]
            num_ents.sort(key=lambda x: x[2], reverse=True)

            if num_ents:
                for f_match in num_ents:
                    matches.append(f_match[0])

            if not all_:
                break

        matches = [match.path if return_type.startswith('file')
                   else match for match in matches]
        return matches if all_ else matches[0] if matches else None

    def get_bvec(self, path, **kwargs):
        """Get bvec file for passed path."""
        result = self.get_nearest(path, extension='bvec', suffix='dwi',
                                  all_=True, **kwargs)
        return listify(result)[0]

    def get_bval(self, path, **kwargs):
        """Get bval file for passed path."""
        result = self.get_nearest(path, suffix='dwi', extension='bval',
                                  all_=True, **kwargs)
        return listify(result)[0]

    def get_fieldmap(self, path, return_list=False):
        """Get fieldmap(s) for specified path."""
        fieldmaps = self._get_fieldmaps(path)

        if return_list:
            return fieldmaps
        else:
            if len(fieldmaps) == 1:
                return fieldmaps[0]
            elif len(fieldmaps) > 1:
                raise ValueError("More than one fieldmap found, but the "
                                 "'return_list' argument was set to False. "
                                 "Either ensure that there is only one "
                                 "fieldmap for this image, or set the "
                                 "'return_list' argument to True and handle "
                                 "the result as a list.")
            else:  # len(fieldmaps) == 0
                return None

    def _get_fieldmaps(self, path):
        sub = self.parse_file_entities(path)['subject']
        fieldmap_set = []
        suffix = '(phase1|phasediff|epi|fieldmap)'
        files = self.get(subject=sub, suffix=suffix, regex_search=True,
                         extension=['nii.gz', 'nii'])
        for file in files:
            metadata = self.get_metadata(file.path)
            if metadata and "IntendedFor" in metadata.keys():
                intended_for = listify(metadata["IntendedFor"])
                if any([path.endswith(_suff) for _suff in intended_for]):
                    cur_fieldmap = {}
                    if file.entities['suffix'] == "phasediff":
                        cur_fieldmap = {"phasediff": file.path,
                                        "magnitude1": file.path.replace(
                                            "phasediff", "magnitude1"),
                                        "suffix": "phasediff"}
                        magnitude2 = file.path.replace(
                            "phasediff", "magnitude2")
                        if os.path.isfile(magnitude2):
                            cur_fieldmap['magnitude2'] = magnitude2
                    elif file.entities['suffix'] == "phase1":
                        cur_fieldmap["phase1"] = file.path
                        cur_fieldmap["magnitude1"] = \
                            file.path.replace("phase1", "magnitude1")
                        cur_fieldmap["phase2"] = \
                            file.path.replace("phase1", "phase2")
                        cur_fieldmap["magnitude2"] = \
                            file.path.replace("phase1", "magnitude2")
                        cur_fieldmap["suffix"] = "phase"
                    elif file.entities['suffix'] == "epi":
                        cur_fieldmap["epi"] = file.path
                        cur_fieldmap["suffix"] = "epi"
                    elif file.entities['suffix'] == "fieldmap":
                        cur_fieldmap["fieldmap"] = file.path
                        cur_fieldmap["magnitude"] = \
                            file.path.replace("fieldmap", "magnitude")
                        cur_fieldmap["suffix"] = "fieldmap"
                    fieldmap_set.append(cur_fieldmap)
        return fieldmap_set

    def get_tr(self, derivatives=False, **filters):

        """Return the scanning repetition time (TR) for one or more runs.

        Parameters
        ----------
        derivatives : bool
            If True, also checks derivatives images.
        filters : dict
            Optional keywords used to constrain the selected runs.
            Can be any arguments valid for a .get call (e.g., BIDS entities
            or JSON sidecar keys).

        Returns
        -------
        float
            A single float.

        Notes
        -----
        Raises an exception if more than one unique TR is found.
        """
        # Constrain search to functional images
        filters.update(suffix='asl', datatype='asl')
        scope = 'all' if derivatives else 'raw'
        images = self.get(extension=['nii', 'nii.gz'], scope=scope,
                          **filters)
        if not images:
            raise ValueError("No functional images that match criteria found.")

        all_trs = set()
        for img in images:
            md = self.get_metadata(img.path)
            all_trs.add(round(float(md['RepetitionTime']), 5))

        if len(all_trs) > 1:
            raise ValueError("Unique TR cannot be found given filters {!r}"
                             .format(filters))
        return all_trs.pop()

    def build_path(self, source, path_patterns=None, strict=False,
                   scope='all', validate=True):
        """Construct a target filename for a file or dictionary of entities.

        Parameters
        ----------
        source : str or :obj:`bids.layout.BIDSFile` or dict
            The source data to use to construct the new file path.
            Must be one of:
            - A BIDSFile object
            - A string giving the path of a BIDSFile contained within the
              current Layout.
            - A dict of entities, with entity names in keys and values in
              values
        path_patterns : list
            Optional path patterns to use to construct
            the new file path. If None, the Layout-defined patterns will
            be used. Entities should be represented by the name
            surrounded by curly braces. Optional portions of the patterns
            should be denoted by square brackets. Entities that require a
            specific value for the pattern to match can pass them inside
            angle brackets. Default values can be assigned by specifying a string
            after the pipe operator. E.g., (e.g., {type<image>|bold} would
            only match the pattern if the entity 'type' was passed and its
            value is "image", otherwise the default value "bold" will be
            used).
                Example: 'sub-{subject}/[var-{name}/]{id}.csv'
                Result: 'sub-01/var-SES/1045.csv'
        strict : bool, optional
            If True, all entities must be matched inside a
            pattern in order to be a valid match. If False, extra entities
            will be ignored so long as all mandatory entities are found.
        scope : str or list, optional
            The scope of the search space. Indicates which
            BIDSLayouts' path patterns to use. See BIDSLayout docstring
            for valid values. By default, uses all available layouts. If
            two or more values are provided, the order determines the
            precedence of path patterns (i.e., earlier layouts will have
            higher precedence).
        validate : bool, optional
            If True, built path must pass BIDS validator. If
            False, no validation is attempted, and an invalid path may be
            returned (e.g., if an entity value contains a hyphen).
        """
        # 'is_file' is a crude check for Path objects
        if isinstance(source, str) or hasattr(source, 'is_file'):
            source = str(source)
            if source not in self.files:
                source = os.path.join(self.root, source)

            source = self.get_file(source)

        if isinstance(source, BIDSFile):
            source = source.entities

        if path_patterns is None:
            layouts = self._get_layouts_in_scope(scope)
            path_patterns = []
            seen_configs = set()
            for l in layouts:
                for c in l.config.values():
                    if c in seen_configs:
                        continue
                    if c.default_path_patterns is not None:
                        path_patterns.extend(c.default_path_patterns)
                    seen_configs.add(c)

        built = build_path(source, path_patterns, strict)
        if built is None:
            raise ValueError(
                "Unable to construct build path with source {}".format(source))
        to_check = os.path.join(os.path.sep, built)

        if not validate or BIDSValidator().is_bids(to_check):
            return built

        raise ValueError("Built path {} is not a valid BIDS filename. Please "
                         "make sure all provided entity values are "
                         "spec-compliant.".format(built))

    def copy_files(self, files=None, path_patterns=None, symbolic_links=True,
                   root=None, conflicts='fail', **kwargs):
        """Copy BIDSFile(s) to new locations.

        The new locations are defined by each BIDSFile's entities and the
        specified `path_patterns`.

        Parameters
        ----------
        files : list
            Optional list of BIDSFile objects to write out. If
            none provided, use files from running a get() query using
            remaining **kwargs.
        path_patterns : str or list
            Write patterns to pass to each file's write_file method.
        symbolic_links : bool
            Whether to copy each file as a symbolic link or a deep copy.
        root : str
            Optional root directory that all patterns are relative
            to. Defaults to current working directory.
        conflicts : str
            Defines the desired action when the output path already exists.
            Must be one of:
                'fail': raises an exception
                'skip' does nothing
                'overwrite': overwrites the existing file
                'append': adds a suffix to each file copy, starting with 1
        kwargs : dict
            Optional key word arguments to pass into a get() query.
        """
        _files = self.get(**kwargs)
        if files:
            _files = list(set(files).intersection(_files))

        for f in _files:
            f.copy(path_patterns, symbolic_link=symbolic_links,
                   root=self.root, conflicts=conflicts)

    def write_contents_to_file(self, entities, path_patterns=None,
                               contents=None, link_to=None,
                               content_mode='text', conflicts='fail',
                               strict=False, validate=True):
        """Write arbitrary data to a file defined by the passed entities and
        path patterns.

        Parameters
        ----------
        entities : dict
            A dictionary of entities, with Entity names in
            keys and values for the desired file in values.
        path_patterns : list
            Optional path patterns to use when building
            the filename. If None, the Layout-defined patterns will be
            used.
        contents : object
            Contents to write to the generate file path.
            Can be any object serializable as text or binary data (as
            defined in the content_mode argument).
        link_to : str
            Optional path with which to create a symbolic link
            to. Used as an alternative to and takes priority over the
            contents argument.
        conflicts : str
            Defines the desired action when the output path already exists.
            Must be one of:
                'fail': raises an exception
                'skip' does nothing
                'overwrite': overwrites the existing file
                'append': adds a suffix to each file copy, starting with 1
        strict : bool
            If True, all entities must be matched inside a
            pattern in order to be a valid match. If False, extra entities
            will be ignored so long as all mandatory entities are found.
        validate : bool
            If True, built path must pass BIDS validator. If
            False, no validation is attempted, and an invalid path may be
            returned (e.g., if an entity value contains a hyphen).
        """
        path = self.build_path(entities, path_patterns, strict,
                               validate=validate)

        if path is None:
            raise ValueError("Cannot construct any valid filename for "
                             "the passed entities given available path "
                             "patterns.")

        write_contents_to_file(path, contents=contents, link_to=link_to,
                               content_mode=content_mode, conflicts=conflicts,
                               root=self.root)


class Query(enum.Enum):
    """Enums for use with BIDSLayout.get()."""
    NONE = 1 # Entity must not be present
    ANY = 2  # Entity must be defined, but with an arbitrary value
