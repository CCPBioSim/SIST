"""Unit tests for sist.ir_finder_cli."""

from __future__ import annotations

import pytest

from sist import ir_finder_cli


class FakeIRFinder:
    """Stand-in for IRFinder that returns a fixed string.

    Attributes:
        temperature: Temperature it was given.
        shape: Shape it was given.
    """

    def __init__(self, *, temperature: float, shape: str) -> None:
        """Record the settings.

        Args:
            temperature: Temperature in kelvin.
            shape: ``"linear"`` or ``"circular"``.
        """
        self.temperature = temperature
        self.shape = shape

    def compute_cruciform_energy_string(self, file_path: str) -> str:
        """Return a fixed inverted repeat string.

        Args:
            file_path: Path to the sequence file.

        Returns:
            The settings and file path, joined by commas.
        """
        return f"{self.temperature},{self.shape},{file_path}"


def test_main_prints_the_ir_string(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(ir_finder_cli, "IRFinder", FakeIRFinder)

    ir_finder_cli.main(["310", "circular", "seq.fa"])

    captured = capsys.readouterr()

    assert captured.out.strip() == "310.0,circular,seq.fa"
    assert "deprecated" in captured.err


def test_main_needs_three_arguments(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as raised:
        ir_finder_cli.main(["310"])

    assert raised.value.code == 1
    assert "usage: IR_finder.pl" in capsys.readouterr().err


def test_main_exits_on_an_error() -> None:
    with pytest.raises(SystemExit) as raised:
        ir_finder_cli.main(["not-a-number", "linear", "seq.fa"])

    assert raised.value.code == 1
