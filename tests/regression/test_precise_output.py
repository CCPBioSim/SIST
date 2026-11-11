"""Regression tests for qsidd's full-precision output (``-P``)."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from conftest import require_success, run_command

from sist.binaries import BinaryLocator

pytestmark = pytest.mark.regression

REGRESSION_DATA = Path(__file__).resolve().parent / "data" / "pbr322.toy.fa"
SETTINGS = ["-T", "310.0", "-s", "0.06", "-i", "0.01", "-t", "12.0"]


@pytest.fixture(scope="session")
def fragment(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The first 2000 bp of the regression sequence, to keep the runs short."""
    lines = REGRESSION_DATA.read_text(encoding="utf-8").splitlines()
    sequence = "".join(line for line in lines if not line.startswith(">"))
    path = tmp_path_factory.mktemp("precise") / "fragment.fa"
    path.write_text(f">fragment\n{sequence[:2000]}\n", encoding="utf-8")

    return path


@pytest.fixture(scope="session")
def binaries(sist_command: list[str]) -> BinaryLocator:
    """Locate the qsidd binaries the rest of the regression suite uses."""
    return BinaryLocator()


def _run_qsidd(binary: Path, fragment: Path, *flags: str) -> str:
    """Run qsidd on the fragment and return its output.

    Args:
        binary: The qsidd binary to run.
        fragment: The sequence file.
        *flags: qsidd flags.

    Returns:
        Standard output.
    """
    process = run_command(
        [str(binary), *flags, *SETTINGS, "-f", fragment.name], cwd=fragment.parent
    )
    require_success(process, description=f"{binary.parent.name} qsidd")

    return process.stdout


def test_precise_output_prints_exact_counts_and_the_log_partition_function(
    binaries: BinaryLocator, fragment: Path
) -> None:
    output = _run_qsidd(binaries.trans_three_binary(), fragment, "-b", "-p", "-r", "-P")

    counts = re.findall(r"Number of \w+-run states = (\S+)", output)

    assert "ln Partition Function = " in output
    assert counts
    assert all(count.isdigit() for count in counts)


def test_precise_competition_output_prints_pairwise_probabilities(
    binaries: BinaryLocator, fragment: Path
) -> None:
    output = _run_qsidd(
        binaries.trans_compete_binary(),
        fragment,
        "-b",
        "-p",
        "-r",
        "-P",
        "-X",
        "1,0,10000,|",
    )

    assert "Joint probability M Z = " in output
    assert "Conditional probability M given Z = " in output
    assert "Conditional probability M given not Z = " in output


def test_output_without_precise_is_unchanged(
    binaries: BinaryLocator, fragment: Path
) -> None:
    output = _run_qsidd(binaries.trans_three_binary(), fragment, "-b", "-p", "-r")

    assert "ln Partition Function" not in output
    assert "Joint probability" not in output
