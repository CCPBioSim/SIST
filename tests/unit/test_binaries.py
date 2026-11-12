"""Unit tests for sist.binaries."""

from __future__ import annotations

from pathlib import Path

import pytest

from sist.binaries import BinaryLocator


def test_binary_uses_the_environment_override(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    binary = tmp_path / "qsidd"
    binary.write_text("", encoding="utf-8")
    monkeypatch.setenv("SIST_TRANS_COMPETE_BIN", str(binary))

    assert BinaryLocator().binary("trans_compete") == binary
    assert BinaryLocator().trans_compete_binary() == binary


def test_binary_reports_a_missing_executable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("SIST_TRANS_THREE_BIN", str(tmp_path / "missing"))

    with pytest.raises(FileNotFoundError, match="SIST_TRANS_THREE_BIN"):
        BinaryLocator().trans_three_binary()


def test_binary_rejects_an_unknown_engine() -> None:
    with pytest.raises(KeyError):
        BinaryLocator().binary("trans_four")
