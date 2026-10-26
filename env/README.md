# Software environment

Record before release:

```r
writeLines(capture.output(sessionInfo()), "env/R_sessionInfo.txt")   # from the machine that ran the R scripts
```
```bash
conda env export --no-builds > env/<name>.yml                        # each conda env on the HPC (pyGenomeTracks, SCENIC+, ldsc, HiCUP)
```
