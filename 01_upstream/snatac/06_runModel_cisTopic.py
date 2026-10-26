import os, sys
sys.path.insert(0, os.path.join(os.environ.get('EHT_REPO', '.'), 'config'))
from paths import DATA_DIR, OUT_DIR, RESOURCE_DIR, TMP_DIR  # noqa: E402
import pickle
import os
import pycisTopic
out_dir = os.path.join(DATA_DIR, "atacseq/cisTopic")


with open(os.path.join(out_dir, "cistopic_obj.pkl"), 'rb') as f:
    cistopic_obj = pickle.load(f)

print(cistopic_obj)




os.environ['MALLET_MEMORY'] = '200G'
from pycisTopic.lda_models import run_cgs_models_mallet
# Configure path Mallet
mallet_path=os.path.join(RESOURCE_DIR, "Mallet-202108/bin/mallet")
# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100],
    n_cpu=12,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=os.path.join(DATA_DIR, "atacseq/tmp"),
    save_path=os.path.join(DATA_DIR, "atacseq/outs"),
    mallet_path=mallet_path,
)

pickle.dump(
    models,
    open(os.path.join(out_dir, "models.pkl"), "wb")
)

