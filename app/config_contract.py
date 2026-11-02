"""Typed, allowlisted configuration contract for local Streamlit runs."""

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable
from urllib.parse import urlparse

ALLOWED_DEPTHS = (5, 10)
ALLOWED_ALIGNERS = ("hisat2",)
REMOTE_SCHEMES = {"http", "https", "gs", "s3", "r2"}


class ConfigValidationError(ValueError):
    """Raised before any subprocess when a run configuration is unsafe."""


@dataclass(frozen=True)
class LocalRunConfig:
    manifest: Path
    output_root: Path
    depth: int = 10
    threads: int = 4
    aligner: str = "hisat2"
    sensitivity_acknowledged: bool = False
    uniqueness_filter: bool = False
    tcrdock: bool = False

    def validate(self, *, require_paths: bool = True) -> "LocalRunConfig":
        errors = []
        if self.depth not in ALLOWED_DEPTHS:
            errors.append(f"depth must be one of {ALLOWED_DEPTHS}")
        if self.depth == 5 and not self.sensitivity_acknowledged:
            errors.append("the 5-read sensitivity preset requires acknowledgement")
        if not 1 <= self.threads <= 8:
            errors.append("threads must be between 1 and 8")
        if self.aligner not in ALLOWED_ALIGNERS:
            errors.append("Phase 0 supports native HISAT2 only")
        if self.uniqueness_filter:
            errors.append("the matched-normal uniqueness filter is locked off")
        if self.tcrdock:
            errors.append("TCRdock execution is outside the local Phase 0 boundary")
        for label, path in (("manifest", self.manifest), ("output_root", self.output_root)):
            if _is_remote(path):
                errors.append(f"{label} must be a local path")
        if require_paths and not self.manifest.is_file():
            errors.append(f"manifest does not exist: {self.manifest}")
        elif require_paths:
            try:
                validate_manifest(self.manifest)
            except ConfigValidationError as exc:
                errors.append(str(exc))
        if errors:
            raise ConfigValidationError("; ".join(errors))
        return self

    def overlay(self) -> Dict[str, Any]:
        self.validate(require_paths=False)
        return {
            "samples_tsv": str(self.manifest),
            "output": {
                "results_dir": str(self.output_root / "results"),
                "logs": str(self.output_root / "logs"),
            },
            "alignment": {
                "aligner": self.aligner,
                "uniqueness_filter": {"enabled": False},
            },
            "filtering": {"min_junction_reads": self.depth},
            "tcrdock": {"enabled": False},
            "streamlit_adapter": {
                "schema_version": 1,
                "sensitivity_acknowledged": self.sensitivity_acknowledged,
            },
        }

    def public_summary(self) -> Dict[str, Any]:
        result = asdict(self)
        result["manifest"] = str(self.manifest)
        result["output_root"] = str(self.output_root)
        return result


def validate_manifest(path: Path) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"patient_id", "sample_type", "fastq1"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ConfigValidationError("manifest is missing columns: " + ", ".join(sorted(missing)))
        rows = list(reader)
    if not rows:
        raise ConfigValidationError("manifest is empty")
    patients = {row["patient_id"].strip() for row in rows if row["patient_id"].strip()}
    if len(patients) != 1:
        raise ConfigValidationError("manifest must contain exactly one patient")
    if not any(row["sample_type"].strip() == "Primary Tumor" for row in rows):
        raise ConfigValidationError("manifest must contain a Primary Tumor sample")
    fastqs = []
    for row in rows:
        for column in ("fastq1", "fastq2"):
            value = row.get(column, "").strip()
            if value:
                fastqs.append(Path(value))
    reject_remote_paths(fastqs)


def _is_remote(path: Path) -> bool:
    text = str(path)
    parsed = urlparse(text)
    return parsed.scheme.lower() in REMOTE_SCHEMES or "://" in text


def reject_remote_paths(paths: Iterable[Path]) -> None:
    remote = [str(path) for path in paths if _is_remote(path)]
    if remote:
        raise ConfigValidationError("remote paths are not allowed: " + ", ".join(remote))
