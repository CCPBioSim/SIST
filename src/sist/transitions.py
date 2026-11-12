"""The structural transitions and algorithms SIST can calculate.

Each transition records what qsidd calls its values. Each algorithm records the
transitions it calculates, and the qsidd binary and flags that run it.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class QsiddNames:
    """What qsidd calls one transition's values in one of its output layouts.

    qsidd has two layouts: one for a single transition (``trans_three``) and
    one for competition between transitions (``trans_compete``).

    Attributes:
        probability: Name of the probability of at least one run.
        expected_runs: Name of the expected number of runs.
        expected_transformed_bp: Name of the expected number of base pairs in
            the transition.
        nucleation_energy: Name of the nucleation energy.
        profile_column: Heading of the per-position probability column.
    """

    probability: str
    expected_runs: str
    expected_transformed_bp: str
    nucleation_energy: str
    profile_column: str


@dataclass(frozen=True)
class Transition:
    """A structural transition SIST can calculate.

    Attributes:
        name: Name used in results, such as ``"z_dna"``.
        letter: Letter qsidd uses for the transition in competition output.
        needs_inverted_repeats: Whether the IR finder must score candidate
            inverted repeats before qsidd runs.
        single: What qsidd calls the values when this is the only transition.
        competition: What qsidd calls the values in competition output.
    """

    name: str
    letter: str
    needs_inverted_repeats: bool
    single: QsiddNames
    competition: QsiddNames

    def names(self, *, competition: bool) -> QsiddNames:
        """Choose the names for one of qsidd's output layouts.

        Args:
            competition: Whether qsidd ran in competition mode.

        Returns:
            What qsidd calls this transition's values in that layout.
        """
        return self.competition if competition else self.single


@dataclass(frozen=True)
class Algorithm:
    """An algorithm the ``sist`` command can run.

    Attributes:
        flag: The ``-a`` flag, such as ``"A"``.
        transitions: Transitions the algorithm calculates, in the order SIST
            reports them.
        engine: qsidd binary that runs it, ``"trans_three"`` or
            ``"trans_compete"``.
        qsidd_flags: Extra flags that pick the transition in the engine.
        supports_nearest_neighbor: Whether ``-n`` is passed on to qsidd.
    """

    flag: str
    transitions: tuple[Transition, ...]
    engine: str
    qsidd_flags: tuple[str, ...]
    supports_nearest_neighbor: bool

    @property
    def transition_names(self) -> tuple[str, ...]:
        """Names of the transitions the algorithm calculates."""
        return tuple(transition.name for transition in self.transitions)

    @property
    def competition(self) -> bool:
        """Whether qsidd runs in competition mode, with its competition layout."""
        return len(self.transitions) > 1

    @property
    def needs_inverted_repeats(self) -> bool:
        """Whether the IR finder must run before qsidd."""
        return any(transition.needs_inverted_repeats for transition in self.transitions)


MELTING = Transition(
    name="melting",
    letter="M",
    needs_inverted_repeats=False,
    single=QsiddNames(
        probability="Transition Probability",
        expected_runs="Average number of runs",
        expected_transformed_bp="Number of melted bases",
        nucleation_energy="Nucleation energy",
        profile_column="P(x)",
    ),
    competition=QsiddNames(
        probability="Prob_M",
        expected_runs="Number of M_runs",
        expected_transformed_bp="N_Melting",
        nucleation_energy="Melting nucleation energy",
        profile_column="P_melt",
    ),
)

Z_DNA = Transition(
    name="z_dna",
    letter="Z",
    needs_inverted_repeats=False,
    single=QsiddNames(
        probability="Transition Probability",
        expected_runs="Average number of runs",
        expected_transformed_bp="Number of Z-DNA bases",
        nucleation_energy="Nucleation energy",
        profile_column="P(x)",
    ),
    competition=QsiddNames(
        probability="Prob_Z",
        expected_runs="Number of Z_runs",
        expected_transformed_bp="N_Z",
        nucleation_energy="Z-DNA nucleation energy",
        profile_column="P_Z",
    ),
)

CRUCIFORM = Transition(
    name="cruciform",
    letter="C",
    needs_inverted_repeats=True,
    single=QsiddNames(
        probability="Transition Probability",
        expected_runs="Average number of runs",
        expected_transformed_bp="Number of cruciform bases",
        nucleation_energy="Nucleation energy",
        profile_column="P(x)",
    ),
    competition=QsiddNames(
        probability="Prob_C",
        expected_runs="Number of C_runs",
        expected_transformed_bp="N_Cruciform",
        nucleation_energy="Cruciform nucleation energy",
        profile_column="P_cruciform",
    ),
)

TRANSITIONS: dict[str, Transition] = {
    transition.name: transition for transition in (MELTING, Z_DNA, CRUCIFORM)
}

ALGORITHMS: dict[str, Algorithm] = {
    algorithm.flag: algorithm
    for algorithm in (
        Algorithm("M", (MELTING,), "trans_three", (), supports_nearest_neighbor=True),
        Algorithm(
            "Z", (Z_DNA,), "trans_three", ("-Z",), supports_nearest_neighbor=False
        ),
        Algorithm(
            "C", (CRUCIFORM,), "trans_three", ("-C",), supports_nearest_neighbor=False
        ),
        Algorithm(
            "A",
            (MELTING, Z_DNA, CRUCIFORM),
            "trans_compete",
            (),
            supports_nearest_neighbor=True,
        ),
    )
}
