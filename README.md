# 3D genome remodeling during in vitro human hematopoietic development is linked to disease trait genetics

Rong Qiu1,^, Khanh B. Trang2,3,^, Carson Shalaby1, James A. Pippin2,3, James Garifallou1, Struan F.A. Grant2,3,4,5,6*, Christopher S. Thom1,6*

1 Division of Neonatology, Children’s Hospital of Philadelphia, Philadelphia, PA, USA

2 Center for Spatial and Functional Genomics, Children's Hospital of Philadelphia, Philadelphia, PA, USA

3 Division of Human Genetics, Children's Hospital of Philadelphia, Philadelphia, PA, USA

4 Department of Genetics, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA, USA

5 Division Endocrinology and Diabetes, Children's Hospital of Philadelphia, Philadelphia, PA, USA

6 Department of Pediatrics, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA, USA

^ These authors contributed equally


*Correspondence

Christopher S. Thom
christopher.thom@pennmedicine.upenn.edu
10-052 Colket Translational Research Building
3501 Civic Center Blvd
Philadelphia, PA 19104

Struan F.A. Grant
grants@pennmedicine.upenn.edu
1102D Abramson Research Center
3615 Civic Center Boulevard
Philadelphia, PA 19104


# Abstract

Three-dimensional genomic architecture provides regulatory scaffolding for cell fate transitions. How chromatin structure underlies hematopoietic development remains poorly understood. Using human induced pluripotent stem cell modeling, we analyzed paired multimodal single nucleus sequencing and chromatin looping to chart chromatin reorganization during in vitro hematopoiesis. Hierarchical reorganization of (sub-)topologically associating domains accompanies endothelial conversion into hematopoietic stem and progenitor cells, revealing structural logic wherein genome topology anticipates formation of cis-regulatory elements that control hematopoietic transcriptional programs. Molecular motor genes are dynamically regulated during hematopoiesis and targeted dynein inhibition promotes hematopoietic progenitor formation, likely through unopposed kinesin activity. Multimodal data integration with genome-wide association study loci reveals enrichment for blood and immune trait variants, connecting in vitro genome remodeling to in vivo human disease genetics. These findings define regulatory architecture of blood cell emergence and establish in vitro hematopoiesis as a tractable platform for interrogating genomic principles underlying human development and disease.

bioRxiv: https://doi.org/10.1101/2025.05.23.655825

## Layout

| Folder | Contents |
|---|---|
| `01_upstream/` | Hi-C tools, versions and parameters (lab-internal pipeline, commands on request); snATAC consensus peaks and pycisTopic (`snatac/`); SCENIC+ configuration and command log (`scenicplus/`); S-LDSC procedure |
| `02_hic_structure/` | HiCRep, TAD / sub-TAD remodeling classes, boundary changes, A/B compartment switching, loop-constrained statistics |
| `03_multiome/` | snATAC per-cell metrics and locus gene-activity statistics |
| `05_gwas/` | S-LDSC plots, variant-to-gene mapping, OpenTargets novelty |
| `06_locus_plots/` | pyGenomeTracks configuration and the ATAC track colour key |
| `config/` | `paths.R` / `paths.py`, the only place paths are set |
| `docs/data_inputs.tsv` | every processed file a script reads, relative to `DATA_DIR` |

Only code that produces a figure panel or table in the final manuscript is included.

## Figure map

| Panel | Script |
|---|---|
| Fig. 2a | `05_gwas/ldsc_plots.R` § *ldsc panel* (bar chart of ATAC peaks and cREs) |
| Figs 2d–f, 3g–h, 4a–b, 4d, 4f–g; Suppl. Figs 10c–d, 10f (track panels) | `06_locus_plots/tracks.ini` (pyGenomeTracks; region set with `--region`) |
| Fig. 3a | `02_hic_structure/hic_structure.R` § *Boundaries aligner analysis* |
| Fig. 3b–c; Suppl. Fig. 9c | `02_hic_structure/hic_structure.R` § *TAD compare* |
| Fig. 3d; Suppl. Fig. 9b | `02_hic_structure/hic_structure.R` § *AB compartment* |
| Fig. 3e | `02_hic_structure/fig3e_plot.R` (plot), `fig3e_stats.R` (exact statistics) |
| Results effect sizes | `02_hic_structure/loop_constrained_effect_sizes.R` |
| Fig. 3g–h colour key | `06_locus_plots/track_color_legend.R` |
| Fig. 6a | `05_gwas/ldsc_plots.R` § *ldsc panel* |
| Fig. 6b; Suppl. Table 12 | `05_gwas/v2g.R` § *output for paper* (inputs from § *re V2g*; OpenTargets novelty from `v2g_opentargets_novelty.R`) |
| Suppl. Fig. 4a | `02_hic_structure/hic_structure.R` § *res 40K HiCRep* |
| Suppl. Fig. 4b | `01_upstream/snatac/` (pycisTopic), `09_Topic_heatmap.R` |
| Suppl. Fig. 5; Suppl. Table 10 | `03_multiome/locus_gene_activity_and_qc.R` |
| Suppl. Fig. 9e–f | `02_hic_structure/plot_loop_constrained.R` [confirm] |
| Fig. 1b–d, 3f | `Scripts` |
| [Fig. 3e input tables `gene_compartment_df.rds`, `region_compartment_df.rds`] | [loop-constrained compartment assignment, to be added] |
| [Fig. 2b–c; Suppl. Fig. 4c; Suppl. Tables 5–7] | [SCENIC+ downstream / motif enrichment scripts, to be added] |
| [Fig. 1e–k, 4c, 4e, 5, 6c–f; Suppl. Figs 1–3, 6–8, 10a–b, 10e, 11, 12; Suppl. Tables 8–9] | [single-cell, bulk RNA-seq and experimental panels, to be added] |
| Suppl. Fig. 9a, 9d | Hi-C pipeline (tools in `01_upstream/README.md`) |

## Running

```bash
export EHT_REPO=$PWD                       # repository root
export EHT_DATA_DIR=/path/to/processed     # layout: docs/data_inputs.tsv
export EHT_OUT_DIR=$PWD/results
Rscript 02_hic_structure/fig3e_stats.R
```

Scripts are the analysis code as run, with only the paths parameterised. They are organised by
section and are meant to be run block by block, not end to end. Several scripts write intermediate
files back into `DATA_DIR` and later scripts read them, keeping the original directory layout.
Locus plots: `cd 06_locus_plots && pyGenomeTracks --tracks tracks.ini --region chr8:17250000-17920000 -o MTUS1.pdf`
(the track files resolve to `../data/` relative to the `.ini`).

## Data availability

Raw and processed Hi-C and Multiome data: GEO [accession]. GWAS summary statistics are public; sources
are listed in `05_gwas/gwas_traits.tsv` (106 traits).

## Software

R 4.x (Seurat, Signac, GenomicRanges, ChIPseeker, ggpubr, pathfindR, data.table); Python 3
(pyGenomeTracks, SCENIC+); HICUP 0.7.4, pairtools 0.3.0, cooler 0.8.11, HiCRep 1.12.2, HiTAD 0.4.2,
TADCompare 1.18.0, Mustache 1.0.1, Fit-Hi-C2 2.0.7, multiHiCcompare 1.8.0, ldsc. Exact versions are in `env/`.

