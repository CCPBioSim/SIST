"""Unit tests for sist.schema_cli."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from sist.schema_cli import main


def _write_result(directory: Path, name: str, result: dict[str, Any]) -> Path:
    """Write a result to a JSON file.

    Args:
        directory: Directory to write the file in.
        name: File name, for example ``"result.json"``.
        result: The result to write.

    Returns:
        Path to the written file.
    """
    result_path = directory / name
    result_path.write_text(json.dumps(result), encoding="utf-8")

    return result_path


def test_main_reports_a_valid_file(
    tmp_path: Path,
    melting_result: dict[str, Any],
    capsys: pytest.CaptureFixture[str],
) -> None:
    result_path = _write_result(tmp_path, "melting.json", melting_result)

    main([str(result_path)])

    assert "melting.json: valid" in capsys.readouterr().out


def test_main_exits_1_for_an_invalid_file(
    tmp_path: Path,
    melting_result: dict[str, Any],
    capsys: pytest.CaptureFixture[str],
) -> None:
    melting_result["calculation"]["algorithm"] = "Z"
    result_path = _write_result(tmp_path, "invalid.json", melting_result)

    with pytest.raises(SystemExit) as raised:
        main([str(result_path)])

    assert raised.value.code == 1
    assert "algorithm Z calculates" in capsys.readouterr().out


def test_main_reports_every_file_before_exiting(
    tmp_path: Path,
    melting_result: dict[str, Any],
    capsys: pytest.CaptureFixture[str],
) -> None:
    missing_path = tmp_path / "missing.json"
    result_path = _write_result(tmp_path, "melting.json", melting_result)

    with pytest.raises(SystemExit):
        main([str(missing_path), str(result_path)])

    output = capsys.readouterr().out

    assert "missing.json: cannot read" in output
    assert "melting.json: valid" in output


def test_main_rejects_nan(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    result_path = tmp_path / "nan.json"
    result_path.write_text('{"value": NaN}', encoding="utf-8")

    with pytest.raises(SystemExit) as raised:
        main([str(result_path)])

    assert raised.value.code == 1
    assert "NaN is not valid JSON" in capsys.readouterr().out


def test_main_exits_2_without_files(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as raised:
        main([])

    assert raised.value.code == 2
    assert "usage: sist-validate" in capsys.readouterr().err


def test_main_exits_1_on_an_unexpected_error(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    def broken_validator() -> None:
        """Fail as an unexpected error would.

        Raises:
            RuntimeError: Always.
        """
        raise RuntimeError("broken")

    monkeypatch.setattr("sist.schema_cli.ResultValidator", broken_validator)

    with pytest.raises(SystemExit) as raised:
        main([str(tmp_path / "result.json")])

    assert raised.value.code == 1
