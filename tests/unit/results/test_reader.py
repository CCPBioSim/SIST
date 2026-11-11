"""Unit tests for sist.results.reader."""

from __future__ import annotations

from results_helpers import OUTPUTS, QSIDD_PATH, read_output

from sist.results.qsidd_output import QsiddOutput
from sist.results.reader import QsiddReader
from sist.transitions import ALGORITHMS


def test_transitions_follow_the_algorithm() -> None:
    for flag, output_name in OUTPUTS.items():
        reader = QsiddReader(ALGORITHMS[flag])

        transitions = reader.transitions(read_output(output_name))

        assert tuple(transitions) == ALGORITHMS[flag].transition_names


def test_transitions_read_the_single_layout(melting_output: QsiddOutput) -> None:
    reader = QsiddReader(ALGORITHMS["M"])

    melting = reader.transitions(melting_output)["melting"]

    assert melting.nucleation_energy_kcal_per_mol == 10.84
    assert melting.probability == melting_output.number("Transition Probability")
    assert len(melting.profile_probability) == 100


def test_transitions_read_the_competition_layout(
    competition_output: QsiddOutput,
) -> None:
    reader = QsiddReader(ALGORITHMS["A"])

    cruciform = reader.transitions(competition_output)["cruciform"]

    assert cruciform.probability == competition_output.number("Prob_C")
    assert cruciform.nucleation_energy_kcal_per_mol == competition_output.number(
        "Cruciform nucleation energy"
    )


def test_free_energy_is_null_where_probability_is_zero(
    melting_output: QsiddOutput,
) -> None:
    reader = QsiddReader(ALGORITHMS["M"])

    melting = reader.transitions(melting_output)["melting"]

    assert melting.free_energy_kcal_per_mol is not None
    assert [value is None for value in melting.free_energy_kcal_per_mol] == [
        value == 0 for value in melting.profile_probability
    ]


def test_free_energy_is_left_out_where_qsidd_does_not_print_it(
    competition_output: QsiddOutput,
) -> None:
    reader = QsiddReader(ALGORITHMS["A"])

    transitions = reader.transitions(competition_output)

    assert all(t.free_energy_kcal_per_mol is None for t in transitions.values())


def test_summary_reads_exact_state_counts(melting_output: QsiddOutput) -> None:
    reader = QsiddReader(ALGORITHMS["M"])

    summary = reader.summary(melting_output)

    assert summary.states_by_number_of_runs == [734, 0, 0, 0]
    assert summary.total_states == 734
    assert summary.pairwise is None


def test_pairwise_is_null_where_undefined() -> None:
    """P(a given b) is undefined when b never happens, so it is null."""
    reader = QsiddReader(ALGORITHMS["A"])

    summary = reader.summary(read_output("competition-no-inverted-repeats.txt"))

    assert summary.pairwise is not None
    assert summary.pairwise["conditional_probability"]["melting"]["cruciform"] is None


def test_given_absent_is_null_when_the_other_transition_is_certain() -> None:
    """qsidd writes 0 for P(a given not b) when b is certain, which is undefined."""
    text = (QSIDD_PATH / "competition.txt").read_text(encoding="utf-8")
    certain = "\n".join(
        "Prob_C = 1" if line.startswith("Prob_C =") else line
        for line in text.splitlines()
    )
    reader = QsiddReader(ALGORITHMS["A"])

    summary = reader.summary(QsiddOutput.parse(certain))

    assert summary.pairwise is not None
    given_absent = summary.pairwise["conditional_probability_given_absent"]
    assert given_absent["melting"]["cruciform"] is None
    assert given_absent["melting"]["z_dna"] is not None


def test_derived_parameters_are_null_where_not_printed() -> None:
    reader = QsiddReader(ALGORITHMS["Z"])

    derived = reader.derived_parameters(read_output("z-dna.txt"))

    assert derived.bat_kcal_per_mol is None
    assert derived.torsional_stiffness_kcal_bp_per_rad2 == 3.6


def test_underflowed() -> None:
    reader = QsiddReader(ALGORITHMS["M"])

    assert reader.underflowed(read_output("underflow.txt"))
    assert not reader.underflowed(read_output("melting.txt"))
