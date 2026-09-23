"""Parity tests for the `sist` Python CLI against the v1.0.0 reference baselines.

These prove the Python implementation is a safe drop-in replacement: it
reuses the exact same tests/reference/v1.0.0/ baselines the existing
regression suite in test_regression.py is already checked against.
"""

from __future__ import annotations

import pytest
from conftest import SistRun
from test_regression import (
    COMPETITION_REFERENCE,
    REFERENCE_DIRECTORY,
    assert_metadata_matches,
    assert_metrics_match,
    parse_competition_output,
    parse_transition_output,
)

pytestmark = pytest.mark.regression


def test_python_competition_command_succeeds(python_competition_run: SistRun) -> None:
    """The Python CLI's competition calculation should complete successfully."""

    process = python_competition_run.process

    assert process.returncode == 0, (
        f"{python_competition_run.name} command failed with exit code "
        f"{process.returncode}\n\nstdout:\n{process.stdout}\n\nstderr:\n{process.stderr}"
    )


def test_python_transition_command_succeeds(python_transition_run: SistRun) -> None:
    """Each Python CLI transition calculation should complete successfully."""

    process = python_transition_run.process

    assert process.returncode == 0, (
        f"{python_transition_run.name} command failed with exit code "
        f"{process.returncode}\n\nstdout:\n{process.stdout}\n\nstderr:\n{process.stderr}"
    )


def test_python_competition_matches_baseline(python_competition_run: SistRun) -> None:
    """The Python CLI's competition output should match the 1.0.0 baseline."""

    expected = parse_competition_output(COMPETITION_REFERENCE)
    actual = parse_competition_output(python_competition_run.output_path)

    assert_metadata_matches(
        name="python competition", expected=expected.metadata, actual=actual.metadata
    )
    assert_metrics_match(
        name="python competition", expected=expected.metrics, actual=actual.metrics
    )

    assert actual.profile.keys() == expected.profile.keys(), (
        "Python competition reported sequence positions changed"
    )

    for position, expected_row in expected.profile.items():
        actual_row = actual.profile[position]
        assert actual_row == expected_row, (
            f"Python competition profile changed at position {position}: "
            f"expected {expected_row!r}, actual {actual_row!r}"
        )


def test_python_transition_matches_baseline(python_transition_run: SistRun) -> None:
    """Each Python CLI transition output should match its 1.0.0 baseline."""

    # Reference files are named after the SIST_TRANSITIONS fixture names.
    name_by_algorithm = {"M": "melting", "Z": "z-dna", "C": "cruciform"}
    reference_output = (
        REFERENCE_DIRECTORY
        / f"{name_by_algorithm[python_transition_run.algorithm]}.txt"
    )

    expected = parse_transition_output(reference_output)
    actual = parse_transition_output(python_transition_run.output_path)

    assert_metadata_matches(
        name=python_transition_run.name,
        expected=expected.metadata,
        actual=actual.metadata,
    )
    assert_metrics_match(
        name=python_transition_run.name,
        expected=expected.metrics,
        actual=actual.metrics,
    )

    assert actual.profile.keys() == expected.profile.keys(), (
        f"{python_transition_run.name} reported sequence positions changed"
    )

    for position, expected_row in expected.profile.items():
        actual_row = actual.profile[position]
        assert actual_row == expected_row, (
            f"{python_transition_run.name} profile changed at position {position}: "
            f"expected {expected_row!r}, actual {actual_row!r}"
        )
