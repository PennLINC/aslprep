"""Functions for working with atlases."""


def get_atlas_names(subset):
    """Get a list of atlases to be used for parcellation and functional connectivity analyses.

    The actual list of files for the atlases is loaded from a different function.

    NOTE: This is a Node function.

    Parameters
    ----------
    subset = {"all", "subcortical", "cortical"}
        Description of the subset of atlases to collect.
        "cortical" atlases are ones with the cortex included, so they can be used for cortical-only
        data (e.g., cortical thickness).
        "subcortical" atlases are ones with *only* subcortical regions.
        "all" combines both lists of atlases.

    Returns
    -------
    :obj:`list` of :obj:`str`
        List of atlases.
    """
    atlases = {
        "cortical": [
            "4S156Parcels",
            "4S256Parcels",
            "4S356Parcels",
            "4S456Parcels",
            "4S556Parcels",
            "4S656Parcels",
            "4S756Parcels",
            "4S856Parcels",
            "4S956Parcels",
            "4S1056Parcels",
            "Glasser",
            "Gordon",
        ],
        "subcortical": [
            "Tian",
            "HCP",
        ],
    }
    atlases["all"] = sorted(list(set(atlases["cortical"] + atlases["subcortical"])))
    return atlases[subset]


def get_atlas_nifti(atlas_name):
    """Select atlas by name from aslprep/data using aslprep.data.load.

    All atlases are in MNI space.

    NOTE: This is a Node function.

    Parameters
    ----------
    atlas_name : {"4S156Parcels", "4S256Parcels", "4S356Parcels", "4S456Parcels", \
                  "4S556Parcels", "4S656Parcels", "4S756Parcels", "4S856Parcels", \
                  "4S956Parcels", "4S1056Parcels", "Glasser", "Gordon", \
                  "Tian", "HCP"}
        The name of the NIFTI atlas to fetch.

    Returns
    -------
    atlas_file : :obj:`str`
        Path to the atlas file.
    atlas_labels_file : :obj:`str`
        Path to the atlas labels file.
    atlas_metadata_file : :obj:`str`
        Path to the atlas metadata file.
    """
    from os.path import isfile, join

    from aslprep.data import load as load_data

    if "4S" in atlas_name or atlas_name in ("Glasser", "Gordon"):
        # 1 mm3 atlases
        atlas_fname = f"tpl-MNI152NLin6Asym_atlas-{atlas_name}_res-01_dseg.nii.gz"
        tsv_fname = f"atlas-{atlas_name}_dseg.tsv"
    else:
        # 2 mm3 atlases
        atlas_fname = f"tpl-MNI152NLin6Asym_atlas-{atlas_name}_res-02_dseg.nii.gz"
        tsv_fname = f"atlas-{atlas_name}_dseg.tsv"

    if "4S" in atlas_name:
        atlas_file = join("/AtlasPack", atlas_fname)
        atlas_labels_file = join("/AtlasPack", tsv_fname)
        atlas_metadata_file = f"/AtlasPack/tpl-MNI152NLin6Asym_atlas-{atlas_name}_dseg.json"
    else:
        atlas_file = load_data.readable(f"atlases/{atlas_fname}").absolute()
        atlas_labels_file = load_data.readable(f"atlases/{tsv_fname}").absolute()
        atlas_metadata_file = load_data.readable(
            f"atlases/tpl-MNI152NLin6Asym_atlas-{atlas_name}_dseg.json",
        ).absolute()

    if not (isfile(atlas_file) and isfile(atlas_labels_file) and isfile(atlas_metadata_file)):
        raise FileNotFoundError(
            f"File(s) DNE:\n\t{atlas_file}\n\t{atlas_labels_file}\n\t{atlas_metadata_file}"
        )

    return atlas_file, atlas_labels_file, atlas_metadata_file
