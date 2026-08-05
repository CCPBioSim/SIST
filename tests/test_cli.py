from __future__ import annotations

import pytest

from conftest import CompetitionRun


pytestmark = pytest.mark.regression


def test_competition_command_succeeds(
    competition_run: CompetitionRun,
) -> None:
    """The documented competition example should complete successfully."""

    process = competition_run.process

    assert process.returncode == 0, (
        f"Competition command failed with exit code "
        f"{process.returncode}\n\n"
        f"stdout:\n{process.stdout}\n\n"
        f"stderr:\n{process.stderr}"
    )


def test_competition_output_is_created(
    competition_run: CompetitionRun,
) -> None:
    """The competition command should create a non-empty output file."""

    output_path = competition_run.output_path

    assert output_path.is_file(), (
        f"Expected output file was not created: {output_path}"
    )

    assert output_path.stat().st_size > 0, (
        f"Output file is empty: {output_path}"
    )


def test_competition_output_contains_expected_sections(
    competition_run: CompetitionRun,
) -> None:
    """The output should contain the expected calculation sections."""

    output = competition_run.output_path.read_text(encoding="utf-8")

    expected_sections = (
        "Sequence Length = 4361",
        "Total Partition Function =",
        "Prob_M =",
        "Prob_Z =",
        "Prob_C =",
        "Position",
        "P_melt",
        "P_Z",
        "P_cruciform",
    )

    for section in expected_sections:
        assert section in output, (
            f"Expected output section was not found: {section!r}"
        )
