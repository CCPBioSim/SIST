"""Unit tests for the example result files in the documentation."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from sist.schema import ResultValidator

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_PATH = REPOSITORY_ROOT / "docs" / "source" / "_static" / "examples" / "results"
EXAMPLE_PATHS = sorted(EXAMPLES_PATH.glob("*.json"))


@pytest.fixture(scope="module")
def validator() -> ResultValidator:
    return ResultValidator()


def test_examples_cover_every_algorithm_and_status() -> None:
    results = [json.loads(path.read_text(encoding="utf-8")) for path in EXAMPLE_PATHS]

    algorithms = {
        result["calculation"]["algorithm"]
        for result in results
        if result["status"] == "success"
    }
    statuses = {result["status"] for result in results}

    assert algorithms == set(ResultValidator.ALGORITHM_TRANSITIONS)
    assert statuses == {"success", "failed"}


@pytest.mark.parametrize("path", EXAMPLE_PATHS, ids=lambda path: path.name)
def test_example_is_valid(validator: ResultValidator, path: Path) -> None:
    result = json.loads(path.read_text(encoding="utf-8"))

    issues = validator.validate(result)

    assert issues == []
    assert result["schema_version"] == "1.0.0"
