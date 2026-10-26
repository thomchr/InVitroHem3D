# Upstream processing

## Hi-C

The Hi-C processing pipeline is internal to the Center for Spatial and Functional Genomics (CHOP).
We report the tools, versions and parameters; the command scripts are available from the
corresponding authors on request.

Two biological replicates per population (adherent D9+7, HPC D9+8), Arima-HiC, hg38.

| Step                                             | Tool                   | Version   | Key parameters                                                                                                                               |
| ------------------------------------------------ | ---------------------- | --------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| Truncation, alignment, filtering, de-duplication | HiCUP (bowtie2)        | 0.7.4     | Arima digest (`hg38_Arima.txt`); lanes processed separately, merged before de-duplication                                                    |
| Read pairs                                       | pairtools / pairix     | 0.3.0 / 0.3.7 | `pairtools parse --assembly hg38` from the HiCUP BAM                                                                                         |
| Contact matrices                                 | cooler                 | 0.8.11    | `cload pairix`; 1, 2, 4, 10, 40 kb, 500 kb, 1 Mb; ICE (`cooler balance`); replicates merged with `cooler merge` for loop calling and display |
| `.hic` files (visualization)                     | juicer_tools `pre`     | 1.22.01   | same resolutions                                                                                                                             |
| Reproducibility                                  | HiCRep (R)             | 1.12.2    | 10 kb, 40 kb, 1 Mb; SCC                                                                                                                      |
| A/B compartments                                 | cooltools `eigs-cis`   | 0.5.1     | 40 kb, GC-content phasing track                                                                                                              |
| Insulation                                       | cooltools `insulation` | 0.5.1     | 10 kb, 500 kb window                                                                                                                         |
| Hierarchical TADs / sub-TADs                     | HiTAD (TADLib)         | 0.4.2     | 10 kb, both replicates as input                                                                                                              |
| TAD comparison                                   | TADCompare (R)         | 1.18.0    | HiTAD domains as `pre_tads`, 10 kb                                                                                                           |
| Loops                                            | Mustache               | 1.0.1     | 1, 2, 4 kb on merged ICE matrices, `-pt 0.1`, 2 Mb distance limit                                                                            |
| Loops                                            | Fit-Hi-C2              | 2.0.7     | 2 kb – 5 Mb, 200 equal-occupancy bins, FDR < 1 × 10⁻⁶                                                                                        |
| Loop merge across resolutions / callers          | lab R script           | —         | —                                                                                                                                            |
| Differential contacts                            | multiHiCcompare / HiCcompare (R) | 1.16.0 / 1.20.0 | 1 kb, fast-lasso normalisation, edgeR GLM; R 4.2.3                                                                                           |
| Loop → OCR / promoter annotation                 | BEDTools `pairToPair`  | 2.30.0    | both anchors; promoters −1,500/+500 bp, GENCODE v40                                                                                          |

## Single-nucleus ATAC: consensus peaks and topic modelling (`snatac/`)

pycisTopic workflow, run in this order: Seurat/Signac ATAC object rebuild (`01`), MACS2 pseudobulk
consensus peaks (`02`), per-cluster QC (`03`), cisTopic object (`04`, `05`), LDA models with Mallet
for 2–100 topics, 500 iterations, α = 50 (by topic), η = 0.1 (`06`), topic binarization and
differentially accessible regions (`07`), plots (`08`, `09`).

## SCENIC+ (`scenicplus/`)

Snakemake pipeline created with `scenicplus init_snakemake`, run with the parameters in `config.yaml`,
which otherwise follow the [SCENIC+ tutorial](https://scenicplus.readthedocs.io/en/latest/tutorials.html)
defaults. The cisTarget motif databases were built from this dataset's consensus peaks
(`create_cisTarget_databases`, v10nr_clust motif collection, 1 kb background padding). Before running,
fill the `${EHT_DATA_DIR}` / `${EHT_RESOURCES}` placeholders, e.g. `envsubst < config.yaml > config.local.yaml`.

Non-default settings: `search_space_upstream` / `search_space_downstream` = `10 500000`
(default `1000 150000`), i.e. region-to-gene links searched up to 500 kb from the gene.

`commands.sh` is the command log for the Multiome analysis folder: cisTarget database build,
SCENIC+ runs, conversion of Seurat objects, and the compartment / TAD boxplots (Fig. 3e). Perturbation
simulations that are not part of the paper were removed.

## S-LDSC

Partitioned heritability followed the procedure in
[mconery/Grant_hFOB_CRISPRi — 2_S-LDSC_Heritability_Enrichment](https://github.com/mconery/Grant_hFOB_CRISPRi/tree/main/2_S-LDSC_Heritability_Enrichment):
cell-type cRE annotations against baseline-LD v2.2, 1000 Genomes Phase 3 EUR LD reference, HapMap3
SNPs. Annotation sets: `narrow`, `atac`, `Prom`, `expand`, `notcRE`, `notcRE_notProm`, per cell cluster.
Plotting and genetic correlation: `05_gwas/ldsc_plots.R`; traits: `05_gwas/gwas_traits.tsv`.

## Multiome alignment and clustering

Cell Ranger ARC [version]; QC, clustering and WNN in Seurat/Signac (`03_multiome/findmarkers.R`,
`snatac/01_remake_ATAC_seurat.R`).
