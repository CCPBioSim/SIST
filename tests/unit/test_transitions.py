"""Unit tests for sist.transitions."""

from __future__ import annotations

import pytest

from sist.schema import ResultValidator
from sist.transitions import ALGORITHMS, TRANSITIONS


def test_transitions_match_the_schema() -> None:
    schema = ResultValidator.load_schema()

    assert set(TRANSITIONS) == set(schema["$defs"]["transition"]["enum"])


def test_algorithms_match_the_schema() -> None:
    schema = ResultValidator.load_schema()

    algorithms = schema["$defs"]["calculation"]["properties"]["algorithm"]["enum"]

    assert set(ALGORITHMS) == set(algorithms)


def test_every_algorithm_uses_known_transitions() -> None:
    for algorithm in ALGORITHMS.values():
        assert set(algorithm.transition_names) <= set(TRANSITIONS)


@pytest.mark.parametrize(
    ("flag", "competition", "needs_inverted_repeats", "engine"),
    [
        ("M", False, False, "trans_three"),
        ("Z", False, False, "trans_three"),
        ("C", False, True, "trans_three"),
        ("A", True, True, "trans_compete"),
    ],
)
def test_algorithm_properties(
    flag: str, competition: bool, needs_inverted_repeats: bool, engine: str
) -> None:
    algorithm = ALGORITHMS[flag]

    assert algorithm.competition is competition
    assert algorithm.needs_inverted_repeats is needs_inverted_repeats
    assert algorithm.engine == engine


def test_names_follow_the_output_layout() -> None:
    melting = TRANSITIONS["melting"]

    single = melting.names(competition=False)
    competition = melting.names(competition=True)

    assert single.profile_column == "P(x)"
    assert competition.profile_column == "P_melt"
    assert competition.probability == "Prob_M"
