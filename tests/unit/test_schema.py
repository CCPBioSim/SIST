"""Unit tests for sist.schema."""

from __future__ import annotations

import copy
from typing import Any

import pytest
from jsonschema import Draft202012Validator

from sist.schema import ResultValidationError, ResultValidator, ValidationIssue


@pytest.fixture(scope="module")
def validator() -> ResultValidator:
    return ResultValidator()


def test_load_schema_finds_the_packaged_schema() -> None:
    schema = ResultValidator.load_schema()

    assert schema["title"] == "SIST result"


def test_load_schema_returns_a_valid_json_schema() -> None:
    schema = ResultValidator.load_schema()

    Draft202012Validator.check_schema(schema)


@pytest.mark.parametrize(
    "fixture_name", ["melting_result", "competition_result", "failed_result"]
)
def test_validate_accepts_a_valid_result(
    validator: ResultValidator,
    request: pytest.FixtureRequest,
    fixture_name: str,
) -> None:
    result = request.getfixturevalue(fixture_name)

    issues = validator.validate(result)

    assert issues == []


def test_validate_points_at_the_bad_value_inside_one_of(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    """A bad value inside a section that can be null is reported where it is."""
    melting_result["transitions"]["melting"]["profile"]["probability"][2] = 1.5

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == [
        "transitions/melting/profile/probability/2"
    ]
    assert "greater than the maximum" in issues[0].message


def test_validate_rejects_results_in_a_failed_run(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["status"] = "failed"

    issues = validator.validate(melting_result)

    assert {"summary", "transitions", "error"} <= {issue.path for issue in issues}


def test_validate_rejects_an_error_in_a_successful_run(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["error"] = {"code": "qsidd_failed", "message": "x", "stage": "qsidd"}

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["error"]


def test_validate_rejects_an_unknown_transition(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    transitions = melting_result["transitions"]
    transitions["b_dna"] = transitions["melting"]

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["transitions"]
    assert "'b_dna' is not one of" in issues[0].message


def test_validate_rejects_major_version_2(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["schema_version"] = "2.0.0"

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["schema_version"]


def test_validate_accepts_a_later_minor_version(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["schema_version"] = "1.7.2"

    issues = validator.validate(melting_result)

    assert issues == []


def test_validate_accepts_unknown_warning_codes(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["warnings"].append(
        {"code": "some_future_code", "message": "x", "source": "sist"}
    )

    issues = validator.validate(melting_result)

    assert issues == []


def test_check_raises_with_every_issue(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["schema_version"] = "2.0.0"
    melting_result["transitions"]["melting"]["profile"]["probability"][2] = 1.5

    with pytest.raises(ResultValidationError) as raised:
        validator.check(melting_result)

    assert [issue.path for issue in raised.value.issues] == [
        "schema_version",
        "transitions/melting/profile/probability/2",
    ]


def test_validation_issue_formats_the_top_level() -> None:
    issue = ValidationIssue("", "'status' is a required property")

    assert str(issue) == "(top level): 'status' is a required property"


def test_validate_rejects_transitions_that_do_not_match_the_algorithm(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["calculation"]["algorithm"] = "Z"

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["calculation/transitions"]
    assert "algorithm Z calculates ['z_dna']" in issues[0].message


def test_validate_rejects_an_extra_transition(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    transitions = melting_result["transitions"]
    transitions["z_dna"] = transitions["melting"]

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["transitions"]
    assert "not in calculation.transitions: ['z_dna']" in issues[0].message


def test_validate_rejects_a_missing_transition(
    validator: ResultValidator, competition_result: dict[str, Any]
) -> None:
    del competition_result["transitions"]["cruciform"]

    issues = validator.validate(competition_result)

    assert [issue.path for issue in issues] == ["transitions"]
    assert "missing transitions: ['cruciform']" in issues[0].message


def test_validate_accepts_free_energy_for_some_transitions(
    validator: ResultValidator, competition_result: dict[str, Any]
) -> None:
    competition_result["transitions"]["z_dna"]["profile"][
        "free_energy_kcal_per_mol"
    ] = [1.0, 2.0, None, 3.0, 4.0]

    issues = validator.validate(competition_result)

    assert issues == []


def test_validate_rejects_a_transition_paired_with_itself(
    validator: ResultValidator, competition_result: dict[str, Any]
) -> None:
    joint_probability = competition_result["summary"]["pairwise"]["joint_probability"]
    joint_probability["melting"] = {**joint_probability["melting"], "melting": 0.1}

    issues = validator.validate(competition_result)

    assert [issue.path for issue in issues] == [
        "summary/pairwise/joint_probability/melting/melting"
    ]


def test_validate_rejects_pairwise_with_one_transition(
    validator: ResultValidator,
    melting_result: dict[str, Any],
    competition_result: dict[str, Any],
) -> None:
    melting_result["summary"]["pairwise"] = copy.deepcopy(
        competition_result["summary"]["pairwise"]
    )

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["summary/pairwise"]


def test_validate_rejects_a_profile_array_of_the_wrong_length(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["transitions"]["melting"]["profile"]["probability"].pop()

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == [
        "transitions/melting/profile/probability"
    ]
    assert "has 4 entries, expected 5" in issues[0].message


def test_validate_rejects_an_analysed_sequence_of_the_wrong_length(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["input"]["analysed_sequence"] += "A"

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["input/analysed_sequence"]


def test_validate_rejects_state_counts_that_do_not_add_up(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    """A count that has lost digits no longer adds up to the total."""
    melting_result["summary"]["state_counts"]["by_number_of_runs"][1] = 4908540

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["summary/state_counts"]


def test_validate_rejects_a_mismatched_ln_partition_function(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["summary"]["ln_partition_function"] += 1.0

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["summary/ln_partition_function"]


def test_validate_rejects_a_zero_partition_function(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["summary"]["partition_function"] = 0.0

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["summary/partition_function"]
    assert "use null" in issues[0].message


def test_validate_accepts_a_null_partition_function(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["summary"]["partition_function"] = None

    issues = validator.validate(melting_result)

    assert issues == []


def test_validate_rejects_a_timestamp_that_is_not_a_date(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["provenance"]["started_at"] = "yesterday"

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["provenance/started_at"]


def test_validate_checks_transitions_in_a_failed_run(
    validator: ResultValidator, failed_result: dict[str, Any]
) -> None:
    failed_result["calculation"]["transitions"] = ["melting"]

    issues = validator.validate(failed_result)

    assert [issue.path for issue in issues] == ["calculation/transitions"]


def test_validate_accepts_a_missing_timestamp(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["provenance"]["started_at"] = None

    issues = validator.validate(melting_result)

    assert issues == []


def test_validate_rejects_pairwise_with_a_transition_not_calculated(
    validator: ResultValidator, competition_result: dict[str, Any]
) -> None:
    competition_result["calculation"]["algorithm"] = "A"
    competition_result["calculation"]["transitions"] = ["melting", "z_dna"]
    del competition_result["transitions"]["cruciform"]

    issues = validator.validate(competition_result)

    assert "uses a transition not in calculation.transitions" in [
        issue.message for issue in issues
    ]


def test_validate_requires_a_sequence_length_for_a_successful_run(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["input"]["sequence_length_bp"] = None

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["input/sequence_length_bp"]
