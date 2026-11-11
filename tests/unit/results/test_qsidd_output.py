"""Unit tests for sist.results.qsidd_output."""

from __future__ import annotations

from results_helpers import read_output

from sist.results.qsidd_output import QsiddOutput


def test_parse_reads_values_that_share_a_line(melting_output: QsiddOutput) -> None:
    assert "BAT" in melting_output.values
    assert "BGC" in melting_output.values


def test_parse_reads_the_scaling_factor(melting_output: QsiddOutput) -> None:
    assert melting_output.number("Scaling factor") == 300.0


def test_parse_reads_warnings(melting_output: QsiddOutput) -> None:
    assert melting_output.warnings == ["sequence length is too short"]


def test_parse_reads_the_profile_table(competition_output: QsiddOutput) -> None:
    assert competition_output.columns == [
        "Position",
        "Base",
        "P_melt",
        "P_Z",
        "P_cruciform",
    ]
    assert len(competition_output.rows) == 100


def test_number_drops_the_unit(melting_output: QsiddOutput) -> None:
    assert melting_output.number("Nucleation energy") == 10.84


def test_number_returns_none_for_nan_and_infinity() -> None:
    output = read_output("underflow.txt")

    assert output.number("Transition Probability") is None
    assert output.number("ln Partition Function") is None


def test_count_reads_an_exact_integer(melting_output: QsiddOutput) -> None:
    assert melting_output.count("Total number of states") == 734


def test_column_reads_one_column(melting_output: QsiddOutput) -> None:
    bases = melting_output.column("Base")

    assert "".join(bases[:5]) == "TTCTC"
