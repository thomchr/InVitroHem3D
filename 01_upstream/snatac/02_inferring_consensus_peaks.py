import os, sys
sys.path.insert(0, os.path.join(os.environ.get('EHT_REPO', '.'), 'config'))
from paths import DATA_DIR, OUT_DIR, RESOURCE_DIR, TMP_DIR  # noqa: E402
import os
import pycisTopic
import pandas as pd
from pybedtools import BedTool


chromsizes = pd.read_table(
     "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
     header = None,
     names = ["Chromosome", "End"]
     )

chromsizes.insert(1, "Start", 0)

macs_out_dir = "macs2_out"


narrow_peak_dict = {}
peak_half_width=250
from pycisTopic.iterative_peak_calling import get_consensus_peaks
path_to_blacklist=os.path.join(RESOURCE_DIR, "TOPMED_hg38_references/hg38_blacklist.bed")

for filename in os.listdir(macs_out_dir):
     # Check if the file ends with ".narrowPeak"
     if filename.endswith('.narrowPeak'):
         # Extract the sample name from the filename (everything before "_peaks.narrowPeak")
         sample_name = filename.split('_peaks.narrowPeak')[0]
         
         # Read in the narrow peak file using pandas or pybedtools
         peaks = pd.read_csv(os.path.join(macs_out_dir, filename), sep='\t', header=0)
         peaks.columns = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'signalValue', 'pvalue', 'qvalue', 'Summit']
         
         # Store the peaks in the dictionary
         narrow_peak_dict[sample_name] = peaks

# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = path_to_blacklist)


consensus_peaks.to_bed(
    path = os.path.join(macs_out_dir, "consensus_regions.bed"),
    keep =True,
    compression = 'infer',
    chain = False)

