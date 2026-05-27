"""Resolve paths from configs/paths.yaml with ${var} interpolation."""

from __future__ import annotations

import os
import re
from functools import lru_cache
from pathlib import Path

import yaml

_VAR_RE = re.compile(r"\$\{([a-zA-Z_][a-zA-Z0-9_]*)\}")


def _interpolate(value: str, flat: dict[str, str]) -> str:
    """Replace ${var} references using a flattened dotless lookup."""
    prev = None
    while prev != value:
        prev = value
        value = _VAR_RE.sub(lambda m: str(flat.get(m.group(1), m.group(0))), value)
    return value


def _flatten(d: dict, prefix: str = "", out: dict | None = None) -> dict:
    out = {} if out is None else out
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            _flatten(v, key, out)
        else:
            out[k] = v
            out[key] = v
    return out


@lru_cache(maxsize=4)
def load_paths(config_path: str | None = None) -> dict:
    """Load paths.yaml with variable interpolation.

    Override location via env var THORACICTCR_PATHS_CONFIG.
    """
    if config_path is None:
        config_path = os.environ.get(
            "THORACICTCR_PATHS_CONFIG",
            str(Path(__file__).resolve().parents[3] / "configs" / "paths.yaml"),
        )
    with open(config_path) as f:
        raw = yaml.safe_load(f)

    flat = _flatten(raw)

    def walk(node):
        if isinstance(node, dict):
            return {k: walk(v) for k, v in node.items()}
        if isinstance(node, str):
            return _interpolate(node, flat)
        return node

    return walk(raw)


def get_path(dotted_key: str, config_path: str | None = None) -> Path:
    """Resolve a dotted key (e.g., 'raw.bam_dir') to a Path."""
    cfg = load_paths(config_path)
    node = cfg
    for part in dotted_key.split("."):
        node = node[part]
    return Path(node)
