dir=${EHT_DATA_DIR}/multiome

mkdir -p $dir/scripts
mkdir -p $dir/log

# Prepare fasta from consensus regions
conda activate create_cistarget_databases
module load BEDTools/2.30.0-GCC-11.3.0

REGION_BED=${EHT_DATA_DIR}/snATACseq/macs2_out/GRanges_clean.bed
GENOME_FASTA="${EHT_RESOURCES}/rgtdata/hg38/genome_hg38.fa"
CHROMSIZES=$dir/hg38.chrom.sizes
SCRIPT_DIR="${EHT_RESOURCES}/create_cisTarget_databases"
${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh ${GENOME_FASTA} ${CHROMSIZES} ${REGION_BED} Hematopoietic_scATAC.hg38.with_1kb_bg_padding.fa 1000 yes

# Create cistarget databases
CBDIR=${EHT_RESOURCES}/aertslab_motif_colleciton/v10nr_clust_public/singletons
FASTA_FILE=$dir/Hematopoietic_scATAC.hg38.with_1kb_bg_padding.fa
ls ${EHT_RESOURCES}/aertslab_motif_colleciton/v10nr_clust_public/singletons > $dir/motifs.txt
sbatch -t 12:00:00 -c 20 --mem 50G -J create_cistarget_motif_databases -o $dir/log/create_cistarget_motif_databases.out $dir/scripts/create_cistarget_motif_databases.sh

# prep RNA-seq
ln -s ${EHT_DATA_DIR}/multiome/HPC_HE_seurat_RNAOnly.h5seurat ./
#R
# anndata rna-seq was throwing error (https://github.com/theislab/scvelo/issues/255)
##fix with adata = ad.read_h5ad("HPC_HE_seurat_RNAOnlynew.h5ad")
#this works (https://github.com/aertslab/scenicplus/issues/381)
python scripts/modify_RNAadata.py

# Run SCENIC+
mkdir -p $dir/outs
mkdir -p $dir/tmp
mkdir -p $dir/scenicplus_outs
mkdir -p $dir/scenicplus_outs500k

conda activate scenicplus
scenicplus init_snakemake --out_dir scplus_pipeline #(default)
scenicplus init_snakemake --out_dir scplus_pipeline2 #(extend the search space up to 500,000)


sbatch -t 24:00:00 -c 10 --mem 200G -J snakemake -o $dir/log/snakemake.out $dir/scripts/run_scenicplus.sh
sbatch -t 24:00:00 -c 10 --mem 200G -J snakemake500k -o $dir/log/snakemake500k.out $dir/scripts/run_scenicplus500k.sh

# downstream 
sbatch  -c 10 --mem 200G -J scenicplus_output -o $dir/log/scenicplus_output.out $dir/scripts/run_python.sh "scenicplus" $dir/scripts scenicplus_output.py

sbatch  -c 10 --mem 200G -J scenicplus_output -o $dir/log/scenicplus_output.out $dir/scripts/run_python.sh "scenicplus" $dir/scripts scenicplus_downstream.py

# plot genes and peaks within ABcompartment, TAD...
sbatch -c 8 --mem 120G -o log/AB_compartment_boxplot.out -J RNA.AB_compartment_boxplot scripts/run_R.sh $dir scripts/ABcompartment_boxplot.R
sbatch -c 8 --mem 120G -o log/DEanalysis_HPC_HE.out -J DEanalysis_HPC_HE scripts/run_R.sh $dir scripts/DEanalysis_HPC_HE.R
# plot TADcomapartment boxplot
sbatch -c 8 --mem 120G -o log/TADcompare_boxplot.out -J HPC_HE.RNA.TADcompare_boxplot scripts/run_R.sh $dir/scripts/TADhier_boxplot.R

# After reviewer, try constrain the peaks by loops
# Get the HE and HPC loops
hicdir=${EHT_DATA_DIR}/hiC/juicer/Hematopoietic_ChrisThom
cut -f1 $hicdir/HE_D9_plus7/merge_loops/mergeRes/HE_D9_plus7.hicLoops.resMerge.FDR_adjusted_withP.anno.gene2OCR.txt | sort -u | sed "s/:/\t/g" > HE_loop.peaks.bed
cut -f1 $hicdir/HPC_Pro_D9_plus8/merge_loops/mergeRes/HPC_Pro_D9_plus8.hicLoops.resMerge.FDR_adjusted_withP.anno.gene2OCR.txt | sort -u | sed "s/:/\t/g" > HPC_loop.peaks.bed

