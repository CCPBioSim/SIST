"""Regression tests for the structured JSON results written by ``sist --json``.

Each algorithm is run on the regression sequence and its result is compared with
the 1.0.0 text references. The references print 6 significant figures while the
JSON holds full precision, so probabilities are compared to that precision.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

import pytest
from conftest import run_command

from sist.schema import ResultValidator

pytestmark = pytest.mark.regression

REGRESSION_PATH = Path(__file__).resolve().parent
REFERENCE_PATH = REGRESSION_PATH / "reference" / "v1.0.0"

JSON_CASES = (
    pytest.param(("melting", "M", "melting.txt", {"melting": "P(x)"}), id="melting"),
    pytest.param(("z-dna", "Z", "z-dna.txt", {"z_dna": "P(x)"}), id="z-dna"),
    pytest.param(
        ("cruciform", "C", "cruciform.txt", {"cruciform": "P(x)"}), id="cruciform"
    ),
    pytest.param(
        (
            "competition",
            "A",
            "competition.rebuilt.txt",
            {"melting": "P_melt", "z_dna": "P_Z", "cruciform": "P_cruciform"},
        ),
        id="competition",
    ),
)

REFERENCE_RTOL = 5e-6


def _read_reference_columns(path: Path) -> dict[str, list[str]]:
    """Read the per-position table from a text reference output.

    Args:
        path: Path to the reference output.

    Returns:
        Each column's values as text, keyed by column name.
    """
    columns: dict[str, list[str]] = {}
    names: list[str] = []

    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("Position"):
            names = line.split("\t")
            columns = {name: [] for name in names}
        elif names:
            for name, value in zip(names, line.split("\t"), strict=False):
                columns[name].append(value)

    return columns


@pytest.fixture(scope="session", params=JSON_CASES)
def json_run(
    request: pytest.FixtureRequest,
    sist_command: list[str],
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[tuple[str, str, str, dict[str, str]], dict[str, Any]]:
    """Run ``sist --json`` once for each algorithm and load the result."""

    case = request.param
    name, algorithm, _, _ = case
    runtime_directory = tmp_path_factory.mktemp(f"sist-json-{name}")
    shutil.copy2(REGRESSION_PATH / "data" / "pbr322.toy.fa", runtime_directory)

    process = run_command(
        [
            *sist_command,
            "-f",
            "pbr322.toy.fa",
            "-a",
            algorithm,
            "--json",
            "result.json",
        ],
        cwd=runtime_directory,
    )

    assert process.returncode == 0, process.stderr

    result_text = (runtime_directory / "result.json").read_text(encoding="utf-8")

    return case, json.loads(result_text)


def test_json_result_is_valid(
    json_run: tuple[tuple[str, str, str, dict[str, str]], dict[str, Any]],
) -> None:
    """Every algorithm should write a result that passes the validator."""

    _, result = json_run

    assert ResultValidator().validate(result) == []
    assert result["status"] == "success"


def test_json_result_matches_the_reference_probabilities(
    json_run: tuple[tuple[str, str, str, dict[str, str]], dict[str, Any]],
) -> None:
    """The per-position probabilities should match the 1.0.0 references."""

    (_, _, reference_name, column_map), result = json_run
    reference = _read_reference_columns(REFERENCE_PATH / reference_name)

    assert result["input"]["analysed_sequence"] == "".join(reference["Base"])

    for transition, column in column_map.items():
        expected = [float(value) for value in reference[column]]
        actual = result["transitions"][transition]["profile"]["probability"]

        assert actual == pytest.approx(expected, rel=REFERENCE_RTOL, abs=0)


def test_json_result_has_exact_state_counts(
    json_run: tuple[tuple[str, str, str, dict[str, str]], dict[str, Any]],
) -> None:
    """State counts should be exact and add up to the total qsidd printed."""

    _, result = json_run
    state_counts = result["summary"]["state_counts"]

    assert sum(state_counts["by_number_of_runs"]) == state_counts["total"]
