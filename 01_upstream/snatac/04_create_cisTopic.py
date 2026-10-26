import os, sys
sys.path.insert(0, os.path.join(os.environ.get('EHT_REPO', '.'), 'config'))
from paths import DATA_DIR, OUT_DIR, RESOURCE_DIR, TMP_DIR  # noqa: E402
import os
import pycisTopic
import os
out_dir = os.path.join(DATA_DIR, "snATACseq/cisTopic")
#        "10_EMP_Ery": os.path.join(DATA_DIR, "atacseq/fragFiles/10_EMP_Ery.bed")
fragments_dict = {
        "HPC": os.path.join(DATA_DIR, "snATACseq/fragFiles/HPC.tsv.gz"),
        "HE": os.path.join(DATA_DIR, "snATACseq/fragFiles/HE.tsv.gz")
        }

from pycisTopic.qc import get_barcodes_passing_qc_for_sample
sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in fragments_dict:
    (
        sample_id_to_barcodes_passing_filters[sample_id],
        sample_id_to_thresholds[sample_id]
    ) = get_barcodes_passing_qc_for_sample(
            sample_id = sample_id,
            pycistopic_qc_output_dir = os.path.join(DATA_DIR, "snATACseq/cisTopic/qc"),
            unique_fragments_threshold = None, # use automatic thresholding
            tss_enrichment_threshold = None, # use automatic thresholding
            frip_threshold = 0,
            use_automatic_thresholds = True,
    )

path_to_regions = os.path.join(os.path.join(DATA_DIR, "snATACseq/macs2_out/GRanges_clean.bed"))
path_to_blacklist = os.path.join(RESOURCE_DIR, "TOPMED_hg38_references/hg38_blacklist.bed")
pycistopic_qc_output_dir = os.path.join(DATA_DIR, "snATACseq/cisTopic/qc")

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl

cistopic_obj_list = []
for sample_id in fragments_dict:
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[sample_id_to_barcodes_passing_filters[sample_id]]

    # Construct the desired barcode format directly
    valid_bc_transformed = [f"{sample_id}_{bc.split('-')[0]}-{bc.split('-')[1]}"
                              for bc in sample_id_to_barcodes_passing_filters[sample_id]]

    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics.reindex(sample_id_to_barcodes_passing_filters[sample_id]), # Ensure metrics align with original barcodes
        valid_bc = sample_id_to_barcodes_passing_filters[sample_id], # Keep original for internal processing
        n_cpu = 10,
        project = sample_id,
        split_pattern = '-'
    )

    # Rename the index of the cistopic_obj *after* creation
    new_index = [f"{sample_id}_{bc.split('-')[0]}-{bc.split('-')[1]}" for bc in cistopic_obj.cell_data.index]
    cistopic_obj.cell_data.index = new_index
    cistopic_obj_list.append(cistopic_obj)
   

import pickle
pickle.dump(
    cistopic_obj_list,
    open(os.path.join(out_dir, "cistopic_obj_list.pkl"), "wb")
)
print(cistopic_obj_list[0])
print(cistopic_obj_list[1])

from pycisTopic.cistopic_class import merge

cistopic_obj = merge(cistopic_obj_list)

# It's possible the merge function might reformat the index.
# If so, we'll need to transform it again after merging.
def transform_merged_index(index):
    parts = index.split('___')
    if len(parts) == 2:
        barcode_sample = parts[0]
        sample = parts[1]
        barcode_parts = barcode_sample.split('-')
        barcode = '-'.join(barcode_parts[:-1])
        leading_sample = barcode_parts[-1]
        return f"{leading_sample}_{barcode}"
    else:
        return index # Handle cases where the index might not have the '___' pattern

new_merged_index = cistopic_obj.cell_data.index.map(transform_merged_index)
cistopic_obj.cell_data.index = new_merged_index

print(cistopic_obj.cell_data.index[:5])


import pickle
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)


