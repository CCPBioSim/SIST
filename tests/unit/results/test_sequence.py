"""Unit tests for sist.results.sequence."""

from __future__ import annotations

from pathlib import Path

import pytest

from sist.results.diagnostics import ErrorCode, ResultError
from sist.results.sequence import SequenceInput


def test_read_skips_the_header_and_uses_upper_case(tmp_path: Path) -> None:
    path = tmp_path / "sequence.fa"
    path.write_text(">example\nacgt\nTTGG\n", encoding="utf-8")

    sequence = SequenceInput.read(path)

    assert sequence.name == "example"
    assert sequence.analysed == "ACGTTTGG"


def test_read_counts_dropped_characters(tmp_path: Path) -> None:
    path = tmp_path / "sequence.fa"
    path.write_text("ACGRYT\n", encoding="utf-8")

    sequence = SequenceInput.read(path)

    assert sequence.analysed == "ACGT"
    assert sequence.dropped == 2


def test_bases_reports_n_as_g_and_finds_its_positions(tmp_path: Path) -> None:
    path = tmp_path / "sequence.fa"
    path.write_text("ACNTN\n", encoding="utf-8")

    sequence = SequenceInput.read(path)

    assert sequence.n_positions == [3, 5]
    assert sequence.bases == "ACGTG"


def test_sha256_matches_the_reported_bases() -> None:
    sequence = SequenceInput(name=None, analysed="TTCTC", dropped=0)

    assert sequence.sha256 == (
        "2f02431f1d80b9dd4c4da090b4aae94406b6129aa7f78d95f6e4f19b3a30c3e2"
    )


@pytest.mark.parametrize(
    ("text", "code"),
    [
        (None, ErrorCode.INPUT_NOT_FOUND),
        (">first\nACGT\n>second\nACGT\n", ErrorCode.MULTIPLE_SEQUENCES),
        (">empty\n", ErrorCode.EMPTY_SEQUENCE),
    ],
)
def test_read_rejects_input_that_cannot_be_analysed(
    tmp_path: Path, text: str | None, code: ErrorCode
) -> None:
    path = tmp_path / "sequence.fa"

    if text is not None:
        path.write_text(text, encoding="utf-8")

    with pytest.raises(ResultError) as raised:
        SequenceInput.read(path)

    assert raised.value.code == code
    assert raised.value.code.stage == "input"
