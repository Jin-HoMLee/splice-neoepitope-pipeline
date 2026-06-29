"""Cloudflare R2 data access for the patient-results notebooks.

Replaces the decommissioned GCS / ``gsutil`` loader (GCP exit, see Issue #854).
Reads each per-result TSV from the project R2 bucket over the S3-compatible API
and caches it locally, so a notebook re-run is offline-fast and reproducible
end-to-end (run all cells from scratch -> identical figures).

Credentials are read from the project-root ``.env`` (gitignored):

    R2_ENDPOINT, R2_BUCKET, R2_ACCESS_KEY_ID, R2_SECRET_ACCESS_KEY

Usage (in a notebook):

    from r2_io import r2_read_tsv_cached
    PREFIX = "results/patient_001"
    report = r2_read_tsv_cached(f"{PREFIX}/reports/report.tsv")

Pass ``refresh=True`` to bypass the cache for a single read. The local cache
lives at ``$R2_CACHE_DIR`` (default ``~/.cache/splice-neoepitope-r2``).
"""

import os
from pathlib import Path

import boto3
import pandas as pd
from botocore.config import Config

_R2_KEYS = ("R2_ENDPOINT", "R2_BUCKET", "R2_ACCESS_KEY_ID", "R2_SECRET_ACCESS_KEY")
_client = None


def _find_project_root(start=None):
    """Walk up from this module until a directory containing ``.env`` is found."""
    here = Path(start or __file__).resolve()
    for parent in [here, *here.parents]:
        if (parent / ".env").exists():
            return parent
    raise FileNotFoundError(
        "Could not locate the project-root .env carrying the R2 credentials."
    )


def _load_env():
    """Populate ``os.environ`` with R2 keys from the project-root ``.env``.

    Existing environment values win (``setdefault``), so an already-sourced
    shell environment is respected. Tolerates ``export KEY=VALUE`` and quotes.
    """
    env_path = _find_project_root() / ".env"
    for raw in env_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        if line.startswith("export "):
            line = line[len("export "):]
        key, _, val = line.partition("=")
        key = key.strip()
        val = val.strip().strip('"').strip("'")
        if key in _R2_KEYS:
            os.environ.setdefault(key, val)
    missing = [k for k in _R2_KEYS if not os.environ.get(k)]
    if missing:
        raise RuntimeError(f"Missing R2 credentials in environment/.env: {missing}")


def _get_client():
    global _client
    if _client is None:
        _load_env()
        _client = boto3.client(
            "s3",
            endpoint_url=os.environ["R2_ENDPOINT"],
            aws_access_key_id=os.environ["R2_ACCESS_KEY_ID"],
            aws_secret_access_key=os.environ["R2_SECRET_ACCESS_KEY"],
            config=Config(signature_version="s3v4", region_name="auto"),
        )
    return _client


def _cache_dir():
    return Path(
        os.environ.get("R2_CACHE_DIR", Path.home() / ".cache" / "splice-neoepitope-r2")
    )


def r2_download_cached(key, *, refresh=False):
    """Download bucket object ``key`` to the local cache; return its local Path.

    Writes to ``<file>.tmp`` then atomically renames, so an interrupted
    download can never poison the cache.
    """
    _load_env()
    bucket = os.environ["R2_BUCKET"]
    local = _cache_dir() / bucket / key
    if local.exists() and not refresh:
        return local
    local.parent.mkdir(parents=True, exist_ok=True)
    tmp = local.with_name(local.name + ".tmp")
    try:
        _get_client().download_file(bucket, key, str(tmp))
        tmp.rename(local)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise
    return local


def r2_read_tsv_cached(key, *, refresh=False, **read_csv_kwargs):
    """Read a TSV object from R2 (cached) into a pandas DataFrame."""
    path = r2_download_cached(key, refresh=refresh)
    return pd.read_csv(path, sep="\t", **read_csv_kwargs)
