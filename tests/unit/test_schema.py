"""Unit tests for sist.schema."""

from __future__ import annotations

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
    melting_result["profile"]["probability"]["melting"][2] = 1.5

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["profile/probability/melting/2"]
    assert "greater than the maximum" in issues[0].message


def test_validate_rejects_results_in_a_failed_run(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["status"] = "failed"

    issues = validator.validate(melting_result)

    assert {"summary", "profile", "error"} <= {issue.path for issue in issues}


def test_validate_rejects_an_error_in_a_successful_run(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    melting_result["error"] = {"code": "qsidd_failed", "message": "x", "stage": "qsidd"}

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["error"]


def test_validate_rejects_an_unknown_transition(
    validator: ResultValidator, melting_result: dict[str, Any]
) -> None:
    probability = melting_result["profile"]["probability"]
    probability["b_dna"] = probability["melting"]

    issues = validator.validate(melting_result)

    assert [issue.path for issue in issues] == ["profile/probability"]
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
    melting_result["profile"]["probability"]["melting"][2] = 1.5

    with pytest.raises(ResultValidationError) as raised:
        validator.check(melting_result)

    assert [issue.path for issue in raised.value.issues] == [
        "schema_version",
        "profile/probability/melting/2",
    ]


def test_validation_issue_formats_the_top_level() -> None:
    issue = ValidationIssue("", "'status' is a required property")

    assert str(issue) == "(top level): 'status' is a required property"
