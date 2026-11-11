"""Saved qsidd outputs and settings shared by the sist.results tests.

The saved qsidd outputs in ``tests/unit/data/qsidd`` were made with ``-b -p -r
-P`` on ``fragment.fa``, the first 20 bp of the regression sequence. The
underflow output was made on the full regression sequence at a superhelical
density of 0.3, with its table cut to 5 rows.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from sist.results.model import RunSettings
from sist.results.qsidd_output import QsiddOutput
from sist.transitions import ALGORITHMS

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
QSIDD_PATH = REPOSITORY_ROOT / "tests" / "unit" / "data" / "qsidd"
FRAGMENT_PATH = QSIDD_PATH / "fragment.fa"
TOY_PATH = REPOSITORY_ROOT / "tests" / "regression" / "data" / "pbr322.toy.fa"
FRAGMENT_IR_STRING = (QSIDD_PATH / "fragment.ir").read_text(encoding="utf-8").strip()

OUTPUTS = {
    "M": "melting.txt",
    "Z": "z-dna.txt",
    "C": "cruciform.txt",
    "A": "competition.txt",
}


def read_output(name: str) -> QsiddOutput:
    """Parse one of the saved qsidd outputs.

    Args:
        name: File name of the saved output, for example ``"melting.txt"``.

    Returns:
        The parsed output.
    """
    return QsiddOutput.parse((QSIDD_PATH / name).read_text(encoding="utf-8"))


def make_settings(flag: str, **overrides: Any) -> RunSettings:
    """Make settings for the 20 bp fragment the saved outputs were run on.

    Args:
        flag: Algorithm flag, such as ``"M"``.
        **overrides: Settings to change from the defaults.

    Returns:
        The settings the saved outputs were made with.
    """
    settings: dict[str, Any] = {
        "file": str(FRAGMENT_PATH),
        "algorithm": ALGORITHMS[flag],
        "temperature": 310.0,
        "superhelical_density": 0.06,
        "salt": 0.01,
        "energy_threshold": 12.0,
        "circular": False,
        "nearest_neighbor": False,
    }
    settings.update(overrides)

    return RunSettings(**settings)
