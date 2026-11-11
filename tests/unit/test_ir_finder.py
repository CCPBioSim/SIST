"""Unit tests for sist.ir_finder."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from sist.ir_finder import IRFinder

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
TOY_PATH = REPOSITORY_ROOT / "tests" / "regression" / "data" / "pbr322.toy.fa"
IRF_REPORT_PATH = (
    REPOSITORY_ROOT
    / "docs"
    / "source"
    / "_static"
    / "examples"
    / "one_line.pbr322.toy.fa.2.10.10.80.10.20.10000.100.1.txt.html"
)


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


@pytest.mark.parametrize("shape", ["linear", "circular"])
def test_compute_cruciform_energy_string_scores_a_saved_irf_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, shape: str
) -> None:
    """Scores pbr322.toy.fa from its saved IRF report, with IRF replaced."""
    monkeypatch.chdir(tmp_path)
    shutil.copy2(TOY_PATH, "pbr322.toy.fa")
    shutil.copy2(IRF_REPORT_PATH, IRF_REPORT_PATH.name)
    monkeypatch.setattr(
        "sist.ir_finder.subprocess.run",
        lambda command, **_: subprocess.CompletedProcess(command, 0),
    )
    finder = IRFinder(temperature=310.0, shape=shape)

    result = finder.compute_cruciform_energy_string("pbr322.toy.fa")

    entries = [entry.split(",")[:2] for entry in result.split("|") if entry]

    assert entries[0] == ["1", "0"]
    assert ["2982", "25"] in entries
    assert ["3089", "48"] in entries
    assert finder.irf_version == "3.05"


def _score_with_reports(
    directory: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    circular_report: str | None,
) -> list[list[str]]:
    """Score pbr322.toy.fa as a circular molecule using saved IRF reports.

    Args:
        directory: Directory to work in.
        monkeypatch: Used to stand in for IRF.
        circular_report: Text of the report for the shifted sequence, or None
            to give the circular pass no report.

    Returns:
        The ``[start, length]`` of every scored entry.
    """
    monkeypatch.chdir(directory)
    shutil.copy2(TOY_PATH, "pbr322.toy.fa")
    shutil.copy2(IRF_REPORT_PATH, IRF_REPORT_PATH.name)

    if circular_report is not None:
        Path(f"circ.{IRF_REPORT_PATH.name}").write_text(
            circular_report, encoding="utf-8"
        )

    monkeypatch.setattr(
        "sist.ir_finder.subprocess.run",
        lambda command, **_: subprocess.CompletedProcess(command, 0),
    )
    result = IRFinder(
        temperature=310.0, shape="circular"
    ).compute_cruciform_energy_string("pbr322.toy.fa")

    return [entry.split(",")[:2] for entry in result.split("|") if entry]


def test_score_skips_loops_longer_than_the_maximum(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    shutil.copy2(TOY_PATH, "pbr322.toy.fa")
    report = IRF_REPORT_PATH.read_text(encoding="utf-8").replace("Loop: 3", "Loop: 101")
    Path(IRF_REPORT_PATH.name).write_text(report, encoding="utf-8")
    monkeypatch.setattr(
        "sist.ir_finder.subprocess.run",
        lambda command, **_: subprocess.CompletedProcess(command, 0),
    )

    result = IRFinder(
        temperature=310.0, shape="linear"
    ).compute_cruciform_energy_string("pbr322.toy.fa")

    starts = [entry.split(",")[:2] for entry in result.split("|") if entry]

    assert ["2982", "25"] not in starts
    assert ["3089", "48"] in starts


def test_circular_pass_maps_repeats_back_to_the_original_positions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Repeats in the shifted sequence move back by half the sequence length."""
    report = IRF_REPORT_PATH.read_text(encoding="utf-8")

    entries = _score_with_reports(tmp_path, monkeypatch, circular_report=report)

    assert ["836", "25"] in entries
    assert ["943", "48"] in entries


def test_circular_pass_skips_repeats_already_found(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A repeat found again across the origin is only scored once."""
    report = IRF_REPORT_PATH.read_text(encoding="utf-8").replace(
        "2982--2992,2996--3006", "837--847,851--861"
    )

    entries = _score_with_reports(tmp_path, monkeypatch, circular_report=report)

    assert entries.count(["2982", "25"]) == 1


@pytest.mark.parametrize(
    ("left_arm", "right_arm", "expected"),
    [
        (["A", "C"], ["*", "-"], (["T", "C", "A", "T"], 8, 2)),
        (["A", "C"], ["*", "*"], (["T", "C", "G"], 9, 1)),
        (["A", "C"], ["*", "A"], (["T", "C", "A"], 9, 1)),
        (["A", "-"], ["*", "G"], (["T", "G", "A", "T"], 8, 3)),
    ],
    ids=["right-gap", "match", "mismatch", "left-gap"],
)
def test_extend_loop_to_minimum(
    left_arm: list[str],
    right_arm: list[str],
    expected: tuple[list[str], int, int],
) -> None:
    """A loop shorter than IRF_MIN_LOOP borrows bases from the ends of the arms."""
    result = IRFinder._extend_loop_to_minimum(["T"], left_arm, right_arm, loop_start=10)

    assert result == expected


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
        (["Inverted Repeats Finder Program written by:", "Version 3.09"], "3.09"),
        (["Sequence: one_line.seq.fa", "Length: 8"], None),
    ],
)
def test_parse_irf_version(lines: list[str], version: str | None) -> None:
    assert IRFinder.parse_irf_version(lines) == version
