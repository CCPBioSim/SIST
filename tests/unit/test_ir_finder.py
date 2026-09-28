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


def test_sequence_length_and_shift_linear(tmp_path: Path) -> None:
    one_line = tmp_path / "one_line.seq.fa"
    one_line.write_text(">one_line.seq.fa\nABCDEFGH", encoding="utf-8")

    length, shifted = IRFinder.sequence_length_and_shift(one_line, circular=False)

    assert length == 8
    assert shifted == ""


def test_sequence_length_and_shift_circular(tmp_path: Path) -> None:
    one_line = tmp_path / "one_line.seq.fa"
    one_line.write_text(">one_line.seq.fa\nABCDEFGH", encoding="utf-8")

    length, shifted = IRFinder.sequence_length_and_shift(one_line, circular=True)

    assert length == 8
    assert shifted == "EFGHABCD"


def test_parse_loop_header_first_example() -> None:
    # From example/one_line.pbr322.toy.fa...txt.html
    line = "    Indices: 2982--2992,2996--3006  Loop: 3  Score: 22"

    loop_length, loop_start, ir_start, ir_end, ir_length = IRFinder.parse_loop_header(
        line
    )

    assert loop_length == 3
    assert loop_start == 2993
    assert ir_start == 2982
    assert ir_end == 3006
    assert ir_length == 25


def test_parse_loop_header_second_example() -> None:
    line = "    Indices: 3089--3110,3115--3136  Loop: 4  Score: 44"

    loop_length, loop_start, ir_start, ir_end, ir_length = IRFinder.parse_loop_header(
        line
    )

    assert loop_length == 4
    assert loop_start == 3111
    assert ir_start == 3089
    assert ir_end == 3136
    assert ir_length == 48


def test_parse_arm_left() -> None:
    line = "       2982 >> AAACCACCGCT >> 2992"

    assert IRFinder.parse_arm(line) == list("AAACCACCGCT")


def test_parse_arm_right() -> None:
    line = "       3006 << *********** << 2996"

    assert IRFinder.parse_arm(line) == list("***********")


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
