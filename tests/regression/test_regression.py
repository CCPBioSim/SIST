"""Regression tests against the SIST 1.0.0 reference outputs.

Every algorithm is run on the regression sequence twice: once with text output
and once with a structured JSON result. The text output must match its
reference exactly, apart from the run time. The JSON result must hold the same
values, compared to the 6 significant figures the references print.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from conftest import SistRun

from sist.results.qsidd_output import QsiddOutput
from sist.transitions import ALGORITHMS

pytestmark = pytest.mark.regression

REFERENCE_DIRECTORY = Path(__file__).resolve().parent / "reference" / "v1.0.0"

REFERENCE_RTOL = 5e-6
"""Relative tolerance for JSON values, which hold more digits than the text."""


def _scientific_lines(text: str) -> list[str]:
    """Split output into lines, leaving out the run time, which always varies.

    Args:
        text: SIST text output.

    Returns:
        The lines that must match the reference.
    """
    return [line for line in text.splitlines() if not line.startswith("Run time")]


def test_text_output_matches_the_reference(text_run: SistRun) -> None:
    reference = REFERENCE_DIRECTORY / text_run.case.reference

    actual = _scientific_lines(text_run.output_path.read_text(encoding="utf-8"))
    expected = _scientific_lines(reference.read_text(encoding="utf-8"))

    assert actual == expected


def test_json_output_matches_the_reference(json_run: SistRun) -> None:
    algorithm = ALGORITHMS[json_run.case.algorithm]
    reference = QsiddOutput.parse(
        (REFERENCE_DIRECTORY / json_run.case.reference).read_text(encoding="utf-8")
    )

    result = json.loads(json_run.output_path.read_text(encoding="utf-8"))

    assert result["input"]["analysed_sequence"] == "".join(reference.column("Base"))
    assert result["summary"]["state_counts"]["total"] == reference.count(
        "Total number of states"
    )

    for transition in algorithm.transitions:
        names = transition.names(competition=algorithm.competition)
        profile = result["transitions"][transition.name]["profile"]
        expected = [float(value) for value in reference.column(names.profile_column)]

        assert profile["probability"] == pytest.approx(
            expected, rel=REFERENCE_RTOL, abs=0
        )

        if "G(x)" in reference.columns:
            free_energy = [
                float(value)
                for value, probability in zip(
                    reference.column("G(x)"), expected, strict=True
                )
                if probability != 0
            ]
            written = [
                value
                for value in profile["free_energy_kcal_per_mol"]
                if value is not None
            ]

            assert written == pytest.approx(free_energy, rel=REFERENCE_RTOL, abs=0)
