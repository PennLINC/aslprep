"""Base Transformation class and associated utilities. """

import re
import warnings
from abc import ABCMeta, abstractmethod
from copy import deepcopy
import itertools
import inspect

import numpy as np
import pandas as pd

from bids.utils import listify
from bids.variables import SparseRunVariable


class Transformation(metaclass=ABCMeta):

    ### Class-level settings ###
    # The following settings govern the way Transformations are applied to the
    # data. The default settings can be overridden within subclasses.

    # List all argument names that specify variables used by Transformation.
    # This is necessary in order to ensure that all and only variables touched
    # by the transformation are cloned before any manipulation occurs.
    # variables in 'variables' are always cloned, so only additional arguments
    # should be specified here.
    _variables_used = ()

    # What data type to pass onto the core _transform() logic. Must be one
    # of 'variable' (the entire BIDSVariable object), 'pandas' (the extracted
    # pandas DF stored in .values), or 'numpy' (just the numpy array inside
    # the .values property of the pandas DF). To minimize overhead and
    # simplify code, it is recommended to avoid using 'variable' if possible.
    _input_type = 'pandas'

    # The data type the internal _transform() method is expected to return.
    # Must be one of 'variable', 'pandas', 'numpy', or 'none'. In the last
    # case, all desired changes must be made in-place within _transform(), as
    # no further changes will be committed.
    _return_type = 'pandas'

    # A tuple indicating which arguments give the names of variables that must
    # all be aligned with one another (i.e., onsets and durations match
    # perfectly) before processing. Defaults to None.
    _align = None

    # Boolean indicating whether the Transformation should be applied to each
    # variable in the input list in turn. When True (default), Transformation
    # is applied once per element in the variable list, with all arguments
    # being passed repeatedly. When False, all data (i.e., variables or their
    # pandas DFs or ndarrays, as specified in _input_type) are passed to the
    # Transformation simultaneously.
    _loopable = True

    # Boolean indicating whether the Transformation can handle groupby
    # operations. When True, a 'groupby' argument is made implicitly available,
    # and if passed, the Transformation will be applied separately to each
    # subset of the data, as defined by the variables named in groupby. When
    # False, the Transformations does not allow grouping, and will raise an
    # exception if groupby is passed. Transformations should set this to False
    # if the groupby argument cannot possibly change the returned result.
    _groupable = True

    # Boolean indicating whether the output argument is mandatory. When False
    # (default), transformations will be applied in-place unless output is set.
    # When True, the user must explicitly specify the output, or an exception
    # is raised.
    _output_required = False

    # An implicit 'dense' argument is always available, and indicates whether
    # or not to operate on dense variables. When True, the arguments listed in
    # _densify control which variables will be densified. Defaults to variables
    # named in the 'variables' argument. Note that if this value is overridden,
    # 'variables' will need to be explicitly included (i.e., the subclass's
    # _densify tuple replaces the base class rather than appending to it).
    _densify = ('variables',)

    # Allow categorical variables in the input arguments? When None (default),
    # any categorical variables encountered as inputs will raise an exception.
    # Otherwise, a tuple giving the names of the arguments whose variables will
    # be passed through as-is even if categorical.
    _allow_categorical = None

    def __new__(cls, collection, variables, *args, **kwargs):
        t = super(Transformation, cls).__new__(cls)
        t._setup(collection, variables, *args, **kwargs)
        return t.transform()

    def _setup(self, collection, variables, *args, **kwargs):
        """Replaces __init__ to set instance attributes because on Python
        >= 3.3, we can't override both new and init. """
        self.collection = collection
        self.variables = listify(variables)
        self.groupby = kwargs.pop('groupby', None)
        self.output = listify(kwargs.pop('output', None))
        self.output_prefix = kwargs.pop('output_prefix', None)
        self.output_suffix = kwargs.pop('output_suffix', None)
        self.dense = kwargs.pop('dense', False)

        # Convert any args to named keyword arguments in order to make sure
        # that operations like densification, alignment, etc. correctly detect
        # all named arguments.
        if args:
            arg_spec = inspect.getfullargspec(self._transform)
            for i, arg_val in enumerate(args):
                # Skip first two argnames--they're always 'self' and
                # 'variables'
                kwargs[arg_spec.args[2 + i]] = arg_val

        self.kwargs = kwargs

        # Expand any detected variable group names or wild cards
        self._expand_variable_groups()
        self._expand_variable_names()

    def _expand_variable_groups(self):
        """ Replace any detected variable groups with the associated lists of
        variable names.
        """
        groups = self.collection.groups
        variables = [groups[v] if v in groups else [v] for v in self.variables]
        self.variables = list(itertools.chain(*variables))

    def _expand_variable_names(self):
        """Filter all available arguments against collection's variables using
        unix-style pattern matching."""
        def _replace_arg_values(values):
            is_iter = isinstance(values, (list, tuple))
            values = listify(values)
            result = []
            # Only try to match strings containing a relevant special character
            for v in values:
                if isinstance(v, str) and re.search('[\*\?\[\]]', v):
                    result.append(self.collection.match_variables(v))
                else:
                    result.append([v])

            result = list(itertools.chain(*result))
            # Don't return a list unless we have to
            if is_iter or len(result) > 1:
                return result
            return result[0]

        # 'variables' is stored separately, so handle it separately
        self.variables = _replace_arg_values(self.variables)

        for k, arg in self.kwargs.items():
            self.kwargs[k] = _replace_arg_values(arg)

    def _clone_variables(self):
        """Deep copy all variables the transformation touches. This prevents us
        from unnecessarily overwriting existing variables. """

        # Always clone the target variables
        self._variables = {v: self.collection[v].clone()
                           for v in self.variables}

        if not self._variables_used:
            return

        # Loop over argument names and clone all variable names in each one
        for var in self._variables_used:
            for v in listify(self.kwargs.get(var, [])):
                # Kludge: we need to allow entity variables to be passed as
                # names even though they don't exist as separate variables
                if (v not in self.collection.variables and
                        v in ['task', 'run', 'session', 'subject']):
                    continue
                self._variables[v] = deepcopy(self.collection[v])

    def _check_categorical_variables(self):
        """Convert categorical variables to dummy-coded indicators. """

        # Collect variable names to pass through
        pass_thru = []
        if self._allow_categorical is not None:
            for arg in self._allow_categorical:
                keys = self.variables if arg == 'variables' \
                    else self.kwargs.get(arg, [])
                pass_thru.extend(listify(keys))
        pass_thru = list(set(pass_thru))

        for name, col in self._variables.items():
            if name not in pass_thru:
                if col.values.values.dtype.kind not in 'bifc':
                    msg = ("The %s transformation does not allow variable '%s'"
                           " to be categorical. Either pass a different "
                           "variable or explicitly convert to a set of binary "
                           "indicators via the 'factor' transformation.")
                    raise ValueError(msg % (self.__class__.__name__, name))

    def _densify_variables(self):

        variables = []

        for var in self._densify:

            if var == 'variables':
                variables.extend(self.variables)
            else:
                variables.extend(listify(self.kwargs.get(var, [])))

        for v in variables:
            var = self._variables[v]
            if isinstance(var, SparseRunVariable):
                sr = self.collection.sampling_rate
                self._variables[v] = var.to_dense(sr)

    def transform(self):

        output_passed = not (self.output is None and self.output_prefix is None
                             and self.output_suffix is None)

        if not output_passed and self._output_required:
            raise ValueError("Transformation '%s' requires output names to be "
                             "provided. Please set at least one of 'output',"
                             "'output_prefix', or 'output_suffix'." %
                             self.__class__.__name__)

        kwargs = self.kwargs

        # Deep copy all variables we expect to touch
        self._clone_variables()

        # Make sure we don't have categorical variables we can't handle
        self._check_categorical_variables()

        # Densify variables if needed
        if self.dense:
            self._densify_variables()

        # Set variables we plan to operate on directly
        variables = [self._variables[c] for c in self.variables]

        # Align variables if needed
        self._align_variables(variables)

        # Pass desired type--variable, DataFrame, or NDArray
        def select_type(col):
            return {'variable': col, 'pandas': col.values,
                    'numpy': col.values.values}[self._input_type]

        data = [select_type(c) for c in variables]

        if not self._loopable:
            variables = [variables]

        for i, col in enumerate(variables):

            # If we still have a list, pass all variables in one block
            if isinstance(col, (list, tuple)):
                result = self._transform(data, **kwargs)
                if self._return_type not in ['none', None]:
                    col = col[0].clone(data=result, name=self.output[0])
            # Otherwise loop over variables individually
            else:
                if self._groupable and self.groupby is not None:
                    result = col.apply(self._transform, groupby=self.groupby,
                                       **kwargs)
                else:
                    result = self._transform(data[i], **kwargs)

            if self._return_type in ['none', None]:
                continue
            elif self._return_type == 'numpy':
                col.values = pd.DataFrame(result)
            elif self._return_type == 'pandas':
                col.values = result
            elif self._return_type == 'variable':
                col = result

            # Overwrite existing variable
            if not output_passed:
                # If multiple variables were returned, add each one separately
                if isinstance(result, (list, tuple)):
                    for r in result:
                        self.collection[r.name] = r
                else:
                    self.collection[col.name] = col

            # Set as a new variable
            else:
                # Either assign new name in order, or re-use existing one
                if self.output is not None:
                    if len(self.variables) == len(self.output) or not \
                            self._loopable:
                        _output = self.output[i]
                    elif len(self.output) == 1:
                        _output = str(self.output) + '_' + col.name
                else:
                    _output = col.name

                # Add prefix and suffix if provided
                if self.output_prefix is not None:
                    _output = self.output_prefix + _output
                if self.output_suffix is not None:
                    _output += self.output_suffix

                col.name = _output
                self.collection[_output] = col

    @abstractmethod
    def _transform(self, **kwargs):
        pass

    def _preprocess(self, col):
        return col

    def _postprocess(self, col):
        return col

    def _align_variables(self, variables, force=True):
        """Checks whether the specified variables have aligned indexes. This
        implies either that all variables are dense, or that all variables are
        sparse and have exactly the same onsets and durations. If variables are
        not aligned and force = True, all variables will be forced to dense
        format in order to ensure alignment.
        """

        if self._align is None or self._align == 'none':
            return

        def _align(variables):
            # If any variable is dense, all variables must be dense
            sparse = [c for c in variables
                      if isinstance(c, SparseRunVariable)]
            if len(sparse) < len(variables):
                if sparse:
                    sparse_names = [s.name for s in sparse]
                    msg = ("Found a mix of dense and sparse variables. May "
                           "cause problems for some transformations.")
                    if force:
                        msg += (" Sparse variables %s were converted to dense "
                                "form to ensure proper alignment." %
                                sparse_names)
                        sr = self.collection.sampling_rate
                        sparse = [s.to_dense(sr) for s in sparse]
                    warnings.warn(msg)
            # If all are sparse, durations, onsets, and index must match
            # perfectly for all
            else:
                def get_col_data(col):
                    return np.c_[col.values.index, col.duration, col.onset]

                def compare_variables(a, b):
                    return len(a) == len(b) and np.allclose(a, b)

                # Compare 1st col with each of the others
                fc = get_col_data(variables[0])
                if not all([compare_variables(fc, get_col_data(c))
                            for c in variables[1:]]):
                    msg = "Misaligned sparse variables found."
                    if force:
                        msg += (" Forcing all sparse variables to dense in "
                                "order to ensure proper alignment.")
                        sr = self.collection.sampling_rate
                        variables = [c.to_dense(sr) for c in variables]
                    warnings.warn(msg)

        align_variables = [listify(self.kwargs[v])
                           for v in listify(self._align) if v in self.kwargs]
        align_variables = list(itertools.chain(*align_variables))
        align_variables = [self.collection[c] for c in align_variables if c]

        if align_variables and self._loopable:
            for c in variables:
                # TODO: should clone all variables in align_variables before
                # alignment to prevent conversion to dense in any given
                # iteration having side effects. This could be an issue if,
                # e.g., some vars in 'variables' are dense and some are sparse.
                _align([c] + align_variables)
        else:
            _align(listify(variables) + align_variables)
