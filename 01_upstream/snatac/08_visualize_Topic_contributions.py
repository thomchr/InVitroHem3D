import os, sys
sys.path.insert(0, os.path.join(os.environ.get('EHT_REPO', '.'), 'config'))
from paths import DATA_DIR, OUT_DIR, RESOURCE_DIR, TMP_DIR  # noqa: E402
import pycisTopic
import os
import pandas as pd

out_dir = os.path.join(DATA_DIR, "snATACseq/cisTopic")

import pickle

with open(os.path.join(out_dir, "cistopic_obj.pkl"), 'rb') as f:
    cistopic_obj = pickle.load(f)





from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)






#    annot_dict[f"pycisTopic_leiden_10_{resolution}"] = {}
#            "cell.type"].value_counts()
#        annot_dict[f"pycisTopic_leiden_10_{resolution}"][cluster] = f"{counts.index[counts.argmax()]}({cluster})"


#    cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'] = [
#        annot_dict[f'pycisTopic_leiden_10_{resolution}'][x] for x in cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'].tolist()

#







my_color_dict = {"cell.type":{"1_Stromal":"#cdc9c9","2_Endo":"#f37f51","3_HE":"#cd5b45","4_HPC":"#f6d706",
    "5_LMPP_Myeloid":"#76eec6", "6_LMPP_Lymphoid":"#66cdaa", "7_EMP":"#f883fa", "8_EMP_Myeloid":"#cd69c9", "9_EMP_Mk":"#8b4789", "10_EMP_Ery":"#cd2890"}} 


cell_type_order = ['1_Stromal','2_Endo','3_HE','4_HPC','5_LMPP_Myeloid','6_LMPP_Lymphoid','7_EMP','8_EMP_Myeloid','9_EMP_Mk','10_EMP_Ery' ]
cistopic_obj.cell_data['cell.type'] = pd.Categorical(cistopic_obj.cell_data['cell.type'], categories=cell_type_order, ordered=True)

cell_topic_heatmap(
    cistopic_obj,
    variables = ['cell.type'] ,
    color_dictionary=my_color_dict,
    scale = False,
    legend_loc_x = 1.1,
    legend_loc_y = -1.2,
    legend_dist_y = -1,
    figsize = (10, 15),
    save = os.path.join(out_dir, "cistopic.contribution.heatmap.pdf")
    )

pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)



