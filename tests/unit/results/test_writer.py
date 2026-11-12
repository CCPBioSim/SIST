"""Unit tests for sist.results.writer."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from results_helpers import make_settings, read_output

from sist.results.builder import ResultBuilder
from sist.results.model import Result
from sist.results.writer import ResultWriter
from sist.schema import ResultValidationError


def _melting_result() -> Result:
    """Build a successful melting result from the saved output.

    Returns:
        The result.
    """
    builder = ResultBuilder(make_settings("M"))
    sequence = builder.read_sequence()

    return builder.success(read_output("melting.txt"), sequence)


def test_to_json_keeps_arrays_of_values_on_one_line() -> None:
    data = _melting_result().to_dict()

    text = ResultWriter.to_json(data)

    assert sum('"probability": [' in line for line in text.splitlines()) == 1
    assert json.loads(text) == data


def test_to_json_rejects_nan() -> None:
    with pytest.raises(ValueError):
        ResultWriter.to_json({"value": float("nan")})


def test_write_writes_a_checked_result(tmp_path: Path) -> None:
    path = tmp_path / "result.json"

    ResultWriter().write(_melting_result(), path)

    assert json.loads(path.read_text(encoding="utf-8"))["status"] == "success"


def test_write_refuses_an_invalid_result(tmp_path: Path) -> None:
    result = _melting_result()
    transitions = result.transitions or {}
    transitions["melting"].profile_probability.append(0.5)

    with pytest.raises(ResultValidationError):
        ResultWriter().write(result, tmp_path / "result.json")

    assert not (tmp_path / "result.json").exists()
