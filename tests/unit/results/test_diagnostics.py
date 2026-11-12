"""Unit tests for sist.results.diagnostics."""

from __future__ import annotations

import pytest

from sist.results.diagnostics import (
    ErrorCode,
    ResultError,
    ResultWarning,
    WarningCode,
    check_sequence,
    check_settings,
    qsidd_warnings,
)
from sist.results.sequence import SequenceInput
from sist.schema import ResultValidator
from sist.transitions import ALGORITHMS


def test_error_stages_are_ones_the_schema_allows() -> None:
    schema = ResultValidator.load_schema()

    allowed = schema["$defs"]["error"]["properties"]["stage"]["enum"]

    assert {code.stage for code in ErrorCode} <= set(allowed)


def test_warning_sources_are_ones_the_schema_allows() -> None:
    schema = ResultValidator.load_schema()

    allowed = schema["$defs"]["warning"]["properties"]["source"]["enum"]

    assert {code.source for code in WarningCode} == set(allowed)


def test_warning_to_dict_includes_the_source_and_details() -> None:
    warning = ResultWarning(WarningCode.SEQUENCE_TOO_SHORT, "Short.", {"bp": 100})

    assert warning.to_dict() == {
        "code": "sequence_too_short",
        "message": "Short.",
        "source": "sist",
        "details": {"bp": 100},
    }


def test_error_to_dict_takes_the_stage_from_the_code() -> None:
    error = ResultError(ErrorCode.QSIDD_FAILED, "Broken.", exit_code=3)

    assert error.to_dict() == {
        "code": "qsidd_failed",
        "message": "Broken.",
        "stage": "qsidd",
        "exit_code": 3,
        "stderr": None,
    }


def test_check_sequence_warns_about_n_bases_dropped_characters_and_length() -> None:
    sequence = SequenceInput(name=None, analysed="ACNTN", dropped=2)

    warnings = check_sequence(sequence)

    assert [warning.code for warning in warnings] == [
        WarningCode.AMBIGUOUS_BASES_REPLACED,
        WarningCode.UNRECOGNISED_CHARACTERS_IGNORED,
        WarningCode.SEQUENCE_TOO_SHORT,
    ]
    assert warnings[0].details == {"count": 2, "positions": [3, 5]}


def test_check_sequence_warns_about_a_long_sequence() -> None:
    sequence = SequenceInput(name=None, analysed="A" * 10001, dropped=0)

    warnings = check_sequence(sequence)

    assert [warning.code for warning in warnings] == [WarningCode.SEQUENCE_TOO_LONG]


@pytest.mark.parametrize(
    ("flag", "nearest_neighbor", "salt", "codes"),
    [
        ("M", True, 0.01, []),
        ("Z", True, 0.01, [WarningCode.NEAREST_NEIGHBOR_IGNORED]),
        ("C", False, 0.1, [WarningCode.CRUCIFORM_SALT_FIXED]),
        ("M", False, 0.1, []),
    ],
)
def test_check_settings(
    flag: str, nearest_neighbor: bool, salt: float, codes: list[WarningCode]
) -> None:
    warnings = check_settings(
        ALGORITHMS[flag], nearest_neighbor=nearest_neighbor, salt=salt
    )

    assert [warning.code for warning in warnings] == codes


def test_qsidd_warnings_skips_the_file_size_length_warning() -> None:
    messages = [
        "threshold is too small, results may be inaccurate",
        "sequence length is too short",
    ]

    warnings = qsidd_warnings(messages)

    assert [warning.code for warning in warnings] == [WarningCode.THRESHOLD_TOO_LOW]
    assert warnings[0].message == "Threshold is too small, results may be inaccurate."
