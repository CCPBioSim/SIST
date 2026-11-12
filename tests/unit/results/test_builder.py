"""Unit tests for sist.results.builder."""

from __future__ import annotations

import pytest
from results_helpers import (
    FRAGMENT_IR_STRING,
    OUTPUTS,
    TOY_PATH,
    make_settings,
    read_output,
)

from sist.results.builder import ResultBuilder
from sist.results.diagnostics import ErrorCode, ResultError, WarningCode
from sist.results.qsidd_output import QsiddOutput
from sist.schema import ResultValidator


@pytest.fixture(scope="module")
def validator() -> ResultValidator:
    return ResultValidator()


@pytest.mark.parametrize(("flag", "output_name"), list(OUTPUTS.items()))
def test_success_builds_a_valid_result(
    validator: ResultValidator, flag: str, output_name: str
) -> None:
    builder = ResultBuilder(make_settings(flag), command=["sist", "-a", flag])
    sequence = builder.read_sequence()

    result = builder.success(
        read_output(output_name),
        sequence,
        ir_string=FRAGMENT_IR_STRING,
        irf_version="3.09",
    )

    assert validator.validate(result.to_dict()) == []
    assert result.status == "success"


def test_success_records_inverted_repeats_only_when_used() -> None:
    builder = ResultBuilder(make_settings("C"))
    sequence = builder.read_sequence()

    result = builder.success(
        read_output("cruciform.txt"), sequence, ir_string=FRAGMENT_IR_STRING
    )

    assert result.inverted_repeats is not None
    assert [repeat.start for repeat in result.inverted_repeats] == [3, 12]


def test_success_warns_when_there_are_no_inverted_repeats() -> None:
    builder = ResultBuilder(make_settings("A"))
    sequence = builder.read_sequence()

    result = builder.success(
        read_output("competition-no-inverted-repeats.txt"),
        sequence,
        ir_string="1,0,10000,|",
    )

    assert WarningCode.NO_INVERTED_REPEATS_FOUND in [w.code for w in result.warnings]


def test_success_replaces_the_qsidd_length_warning_with_its_own() -> None:
    """Only SIST's base pair length warning is recorded, not qsidd's."""
    builder = ResultBuilder(make_settings("M"))
    sequence = builder.read_sequence()

    result = builder.success(read_output("melting.txt"), sequence)

    assert [w.code for w in result.warnings] == [WarningCode.SEQUENCE_TOO_SHORT]


def test_underflow_becomes_a_failed_result_that_keeps_what_qsidd_printed(
    validator: ResultValidator,
) -> None:
    builder = ResultBuilder(
        make_settings("M", file=str(TOY_PATH), superhelical_density=0.3)
    )
    sequence = builder.read_sequence()
    output = read_output("underflow.txt")

    with pytest.raises(ResultError) as raised:
        builder.success(output, sequence)

    result = builder.failure(raised.value, sequence=sequence, output=output)

    assert raised.value.code == ErrorCode.PARTITION_FUNCTION_UNDERFLOW
    assert validator.validate(result.to_dict()) == []
    assert result.derived_parameters is not None
    assert [w.code for w in result.warnings] == [
        WarningCode.SUPERHELICAL_DENSITY_OUT_OF_RANGE
    ]


def test_failure_before_the_sequence_is_read_is_valid(
    validator: ResultValidator,
) -> None:
    builder = ResultBuilder(make_settings("A", file="missing.fa"))
    error = ResultError(ErrorCode.INPUT_NOT_FOUND, "Missing.")

    result = builder.failure(error)

    assert validator.validate(result.to_dict()) == []
    assert result.status == "failed"
    assert result.input.sequence_length_bp is None


def test_failure_leaves_out_derived_parameters_it_cannot_read() -> None:
    builder = ResultBuilder(make_settings("M"))
    error = ResultError(ErrorCode.INTERNAL_ERROR, "Unreadable.")

    result = builder.failure(error, output=QsiddOutput())

    assert result.derived_parameters is None
