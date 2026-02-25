"""Functions for generating boilerplate test."""

from aslprep.utils.misc import list_to_str


def describe_atlases(atlases):
    """Build a text description of the atlases that will be used."""
    atlas_descriptions = {
        'Glasser': 'the Glasser atlas [@Glasser_2016]',
        'Gordon': 'the Gordon atlas [@Gordon_2014]',
        'Tian': 'the Tian subcortical atlas [@tian2020topographic]',
        'HCP': 'the HCP CIFTI subcortical atlas [@glasser2013minimal]',
        'MIDB': (
            'the MIDB precision brain atlas derived from ABCD data and thresholded at 75% '
            'probability [@hermosillo2024precision]'
        ),
        'MyersLabonte': (
            'the Myers-Labonte infant atlas thresholded at 50% probability [@myers2023functional]'
        ),
    }

    atlas_strings = []
    described_atlases = []
    atlases_4s = [atlas for atlas in atlases if str(atlas).startswith('4S')]
    described_atlases += atlases_4s
    if atlases_4s:
        parcels = [int(str(atlas[2:-7])) for atlas in atlases_4s]
        s = (
            'the Schaefer Supplemented with Subcortical Structures (4S) atlas '
            '[@Schaefer_2017;@pauli2018high;@king2019functional;@najdenovska2018vivo;'
            '@glasser2013minimal] '
            f'at {len(atlases_4s)} different resolutions ({list_to_str(parcels)} parcels)'
        )
        atlas_strings.append(s)

    for k, v in atlas_descriptions.items():
        if k in atlases:
            atlas_strings.append(v)
            described_atlases.append(k)

    undescribed_atlases = [atlas for atlas in atlases if atlas not in described_atlases]
    for atlas in undescribed_atlases:
        atlas_strings.append(f'the {atlas} atlas')

    return list_to_str(atlas_strings)
