"""Unit tests for sist.results.model."""

from __future__ import annotations

from results_helpers import make_settings

from sist.results.model import InputDescription, InvertedRepeat
from sist.results.sequence import SequenceInput


def test_parameters_record_the_superhelical_density_as_used() -> None:
    parameters = make_settings("M").parameters()

    assert parameters["superhelical_density"] == -0.06
    assert parameters["melting_energetics"] == "copolymeric"
    assert "cruciform_salt_concentration_m" not in parameters


def test_parameters_record_the_fixed_cruciform_salt() -> None:
    parameters = make_settings("C", salt=0.1).parameters()

    assert parameters["salt_concentration_m"] == 0.1
    assert parameters["cruciform_salt_concentration_m"] == 0.01
    assert parameters["melting_energetics"] is None


def test_parameters_record_nearest_neighbor_energetics() -> None:
    parameters = make_settings("A", nearest_neighbor=True).parameters()

    assert parameters["melting_energetics"] == "nearest_neighbor"


def test_input_description_of_an_unread_sequence_leaves_out_the_bases() -> None:
    description = InputDescription.describe(make_settings("M"), None)

    section = description.to_dict()

    assert section["sequence_length_bp"] is None
    assert "analysed_sequence" not in section


def test_input_description_reports_the_analysed_bases() -> None:
    sequence = SequenceInput(name="example", analysed="ACNT", dropped=0)

    section = InputDescription.describe(make_settings("M"), sequence).to_dict()

    assert section["analysed_sequence"] == "ACGT"
    assert section["sequence_length_bp"] == 4


def test_parse_all_skips_the_placeholder_entry() -> None:
    inverted_repeats = InvertedRepeat.parse_all("1,0,10000,|20,25,1.5,|60,12,2.5,|")

    assert inverted_repeats == [
        InvertedRepeat(20, 25, 1.5),
        InvertedRepeat(60, 12, 2.5),
    ]


def test_parse_all_of_an_empty_string() -> None:
    assert InvertedRepeat.parse_all("") == []
