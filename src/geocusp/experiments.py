from __future__ import annotations

import json
import platform
import subprocess
from pathlib import Path

import numpy as np

from .backends import detect


def env_report() -> dict:
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        commit = "not-a-git-checkout"
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "git_commit": commit,
        "backends": detect(),
    }


def write_env_report(path: str = "results/environment.json") -> dict:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    out = env_report()
    Path(path).write_text(json.dumps(out, indent=2))
    return out


def save_json(data: dict, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, default=str))
