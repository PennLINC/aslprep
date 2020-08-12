from collections import namedtuple
from aslprep.niworkflows.utils.spaces import Reference, SpatialReferences
from aslprep.workflows.base import init_single_subject_wf
BIDSLayout = namedtuple('BIDSLayout', ('root'))
wf = init_single_subject_wf('01')