import os, sys
sys.path.insert(0, os.path.join(os.environ.get('EHT_REPO', '.'), 'config'))
from paths import DATA_DIR, OUT_DIR, RESOURCE_DIR, TMP_DIR  # noqa: E402
import pickle
import os
import pycisTopic
out_dir = os.path.join(DATA_DIR, "atacseq/cisTopic")

with open(os.path.join(out_dir, "cistopic_obj_list.pkl"), 'rb') as f:
    cistopic_obj_list = pickle.load(f)

from pycisTopic.cistopic_class import merge

cistopic_obj = merge(cistopic_obj_list)

print(cistopic_obj)

pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)
