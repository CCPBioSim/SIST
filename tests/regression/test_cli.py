from __future__ import annotations

import pytest

from conftest import SistRun


pytestmark = pytest.mark.regression


def assert_command_succeeds(sist_run: SistRun) -> None:
    """Assert that a SIST calculation completed successfully."""

    process = sist_run.process

    assert process.returncode == 0, (
        f"{sist_run.name} command failed with exit code "
        f"{process.returncode}\n\n"
        f"stdout:\n{process.stdout}\n\n"
        f"stderr:\n{process.stderr}"
    )


def assert_output_is_created(sist_run: SistRun) -> None:
    """Assert that a SIST calculation created a non-empty output file."""

    output_path = sist_run.output_path

    assert output_path.is_file(), (
        f"Expected output file was not created: {output_path}"
    )

    assert output_path.stat().st_size > 0, (
        f"Output file is empty: {output_path}"
    )


def test_competition_command_succeeds(
    competition_run: SistRun,
) -> None:
    """The competition calculation should complete successfully."""

    assert_command_succeeds(competition_run)


def test_competition_output_is_created(
    competition_run: SistRun,
) -> None:
    """The competition calculation should create a non-empty output file."""

    assert_output_is_created(competition_run)


def test_competition_output_contains_expected_sections(
    competition_run: SistRun,
) -> None:
    """Competition output should contain the expected calculation sections."""

    output = competition_run.output_path.read_text(encoding="utf-8")

    expected_sections = (
        "Sequence Length =",
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


def test_transition_command_succeeds(
    transition_run: SistRun,
) -> None:
    """Each individual transition calculation should complete successfully."""

    assert_command_succeeds(transition_run)


def test_transition_output_is_created(
    transition_run: SistRun,
) -> None:
    """Each transition calculation should create a non-empty output file."""

    assert_output_is_created(transition_run)


def test_transition_output_contains_expected_sections(
    transition_run: SistRun,
) -> None:
    """Transition output should contain the expected calculation sections."""

    output = transition_run.output_path.read_text(encoding="utf-8")

    expected_sections = (
        "Sequence Length =",
        "Total Partition Function =",
        "Transition Probability =",
        "Position",
        "P(x)",
    )

    for section in expected_sections:
        assert section in output, (
            f"{transition_run.name}: expected output section "
            f"was not found: {section!r}"
        )
