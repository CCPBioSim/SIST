"""Unit tests for sist.ir_finder."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from sist.ir_finder import IRFinder


def test_convert_to_one_line_keeps_existing_header(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    Path("seq.fa").write_text(">seq.fa\nACGT\nTTTT\n", encoding="utf-8")

    output_path = IRFinder.convert_to_one_line(Path("seq.fa"))

    assert output_path == Path("one_line.seq.fa")
    assert output_path.read_text(encoding="utf-8") == ">seq.fa\nACGTTTTT"


def test_convert_to_one_line_discards_first_line_without_header(
    tmp_path: Path, monkeypatch
) -> None:
    """A headerless first line is dropped entirely, not kept as sequence data."""
    monkeypatch.chdir(tmp_path)
    Path("seq.fa").write_text("DROPPED_LINE\nACGT\nTTTT\n", encoding="utf-8")

    output_path = IRFinder.convert_to_one_line(Path("seq.fa"))

    assert output_path.read_text(encoding="utf-8") == ">one_line.seq.fa\nACGTTTTT"


def test_convert_to_one_line_empty_file(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    Path("seq.fa").write_text("", encoding="utf-8")

    output_path = IRFinder.convert_to_one_line(Path("seq.fa"))

    assert output_path.read_text(encoding="utf-8") == ""


def test_read_sequence(tmp_path: Path) -> None:
    one_line = tmp_path / "one_line.seq.fa"
    one_line.write_text(">one_line.seq.fa\nABCDEFGH", encoding="utf-8")

    assert IRFinder.read_sequence(one_line) == "ABCDEFGH"


@pytest.mark.parametrize(("circular", "shifted"), [(False, ""), (True, "EFGHABCD")])
def test_sequence_length_and_shift(
    tmp_path: Path, circular: bool, shifted: str
) -> None:
    one_line = tmp_path / "one_line.seq.fa"
    one_line.write_text(">one_line.seq.fa\nABCDEFGH", encoding="utf-8")

    length, actual = IRFinder.sequence_length_and_shift(one_line, circular=circular)

    assert length == 8
    assert actual == shifted


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        (
            "    Indices: 2982--2992,2996--3006  Loop: 3  Score: 22",
            (3, 2993, 2982, 3006, 25),
        ),
        (
            "    Indices: 3089--3110,3115--3136  Loop: 4  Score: 44",
            (4, 3111, 3089, 3136, 48),
        ),
    ],
)
def test_parse_loop_header(line: str, expected: tuple[int, ...]) -> None:
    """The examples come from the IRF report for pbr322.toy.fa in the docs."""
    assert IRFinder.parse_loop_header(line) == expected


@pytest.mark.parametrize(
    ("line", "arm"),
    [
        ("       2982 >> AAACCACCGCT >> 2992", "AAACCACCGCT"),
        ("       3006 << *********** << 2996", "***********"),
    ],
)
def test_parse_arm(line: str, arm: str) -> None:
    assert IRFinder.parse_arm(line) == list(arm)


@pytest.mark.skipif(shutil.which("irf") is None, reason="requires the irf executable")
def test_compute_cruciform_energy_string_smoke(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    Path("toy.fa").write_text(
        ">toy.fa\nAAACCACCGCTTTTTTTTTTTTTTTTTTTTTTTTTGGCGGTGGTTT\n",
        encoding="utf-8",
    )

    finder = IRFinder(temperature=310.0, shape="linear")
    result = finder.compute_cruciform_energy_string("toy.fa")

    assert result.startswith("1,0,10000,|")


@pytest.mark.parametrize(
    ("lines", "version"),
    [
        (["Inverted Repeats Finder Program written by:", "Version 3.05"], "3.05"),
        (["Sequence: one_line.seq.fa", "Length: 8"], None),
    ],
)
def test_parse_irf_version(lines: list[str], version: str | None) -> None:
    assert IRFinder.parse_irf_version(lines) == version
