"""Shared paths; see config/paths.R."""
import os

REPO_DIR = os.path.abspath(os.environ.get("EHT_REPO", "."))
DATA_DIR = os.path.abspath(os.environ.get("EHT_DATA_DIR", os.path.join(REPO_DIR, "data")))
OUT_DIR = os.path.abspath(os.environ.get("EHT_OUT_DIR", os.path.join(REPO_DIR, "results")))
os.makedirs(OUT_DIR, exist_ok=True)
RESOURCE_DIR = os.path.abspath(os.environ.get("EHT_RESOURCES", os.path.join(REPO_DIR, "resources")))
TMP_DIR = os.environ.get("TMPDIR", "/tmp")
