# Shared paths. Set EHT_DATA_DIR to the directory holding the processed data
# (layout in README.md, "Data layout"), and EHT_OUT_DIR for figures and tables.
REPO_DIR <- normalizePath(Sys.getenv("EHT_REPO", "."), mustWork = FALSE)
DATA_DIR <- normalizePath(Sys.getenv("EHT_DATA_DIR", file.path(REPO_DIR, "data")), mustWork = FALSE)
OUT_DIR  <- normalizePath(Sys.getenv("EHT_OUT_DIR",  file.path(REPO_DIR, "results")), mustWork = FALSE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
RESOURCE_DIR <- normalizePath(Sys.getenv("EHT_RESOURCES", file.path(REPO_DIR, "resources")), mustWork = FALSE)  # genome, blacklist, motif collections
