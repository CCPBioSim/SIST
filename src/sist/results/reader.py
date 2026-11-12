"""Turn qsidd's output into result values, for any algorithm.

qsidd prints each value under one name for a single transition and another for
a competition. The reader looks the names up in :mod:`sist.transitions` for the
algorithm it was given.
"""

from __future__ import annotations

from sist.results.model import (
    DerivedParameters,
    Diagnostics,
    PairMatrices,
    Summary,
    TransitionResult,
)
from sist.results.qsidd_output import QsiddOutput
from sist.transitions import Algorithm, Transition


class QsiddReader:
    """Read result values from the full-precision output of one qsidd run.

    Attributes:
        FREE_ENERGY_COLUMN: Heading of qsidd's G(x) column, which it prints for
            a single melting transition.
        algorithm: The algorithm qsidd ran.
    """

    FREE_ENERGY_COLUMN = "G(x)"

    def __init__(self, algorithm: Algorithm) -> None:
        """Set up a reader for one algorithm's output.

        Args:
            algorithm: The algorithm qsidd ran.
        """
        self.algorithm = algorithm

    def underflowed(self, output: QsiddOutput) -> bool:
        """Check whether the partition function underflowed to zero.

        qsidd reports this by printing a log of the partition function that is
        not a finite number.

        Args:
            output: The parsed qsidd output.

        Returns:
            True if the partition function underflowed.
        """
        return output.number("ln Partition Function") is None

    def derived_parameters(self, output: QsiddOutput) -> DerivedParameters:
        """Read the model values qsidd works out from the settings.

        Args:
            output: The parsed qsidd output.

        Returns:
            The derived parameters.

        Raises:
            KeyError: If qsidd did not print them.
        """
        stiffness_key = "K" if "K" in output.values else "C"

        return DerivedParameters(
            linking_difference_turns=output.number("linking difference"),
            torsional_stiffness_kcal_bp_per_rad2=output.number(stiffness_key),
            bat_kcal_per_mol=self._optional_number(output, "BAT"),
            bgc_kcal_per_mol=self._optional_number(output, "BGC"),
        )

    def summary(self, output: QsiddOutput) -> Summary:
        """Read the whole-molecule results.

        Args:
            output: The parsed qsidd output.

        Returns:
            The summary.
        """
        return Summary(
            ln_partition_function=output.number("ln Partition Function"),
            partition_function=output.number("Total Partition Function") or None,
            average_number_of_runs=output.number("Average number of runs"),
            states_by_number_of_runs=[
                output.count(f"Number of {runs}-run states")
                for runs in ("one", "two", "three", "four")
            ],
            total_states=output.count("Total number of states"),
            pairwise=self._pairwise(output) if self.algorithm.competition else None,
        )

    def transitions(self, output: QsiddOutput) -> dict[str, TransitionResult]:
        """Read the results for each transition.

        Args:
            output: The parsed qsidd output.

        Returns:
            The results, keyed by transition name.
        """
        return {
            transition.name: self._transition(output, transition)
            for transition in self.algorithm.transitions
        }

    def diagnostics(self, output: QsiddOutput) -> Diagnostics:
        """Read the values for checking a run.

        Args:
            output: The parsed qsidd output.

        Returns:
            The diagnostics.
        """
        return Diagnostics(
            scaling_factor=output.number("Scaling factor"),
            min_energy_kcal_per_mol=output.number("MinE"),
            max_energy_kcal_per_mol=output.number("Max_E"),
        )

    def _transition(
        self, output: QsiddOutput, transition: Transition
    ) -> TransitionResult:
        """Read the results for one transition.

        G(x) is included when qsidd printed it, with null where the
        transition's probability is 0.

        Args:
            output: The parsed qsidd output.
            transition: The transition to read.

        Returns:
            The transition's results.
        """
        names = transition.names(competition=self.algorithm.competition)
        probability = [float(value) for value in output.column(names.profile_column)]
        free_energy = None

        if self.FREE_ENERGY_COLUMN in output.columns:
            free_energy = [
                None if position_probability == 0 else float(value)
                for position_probability, value in zip(
                    probability, output.column(self.FREE_ENERGY_COLUMN), strict=True
                )
            ]

        return TransitionResult(
            nucleation_energy_kcal_per_mol=output.number(names.nucleation_energy),
            probability=output.number(names.probability),
            expected_runs=output.number(names.expected_runs),
            expected_transformed_bp=output.number(names.expected_transformed_bp),
            profile_probability=probability,
            free_energy_kcal_per_mol=free_energy,
        )

    def _pairwise(self, output: QsiddOutput) -> PairMatrices:
        """Read the probabilities for pairs of transitions.

        P(a given not b) is null when b is certain. Any value qsidd prints as
        ``nan`` is null.

        Args:
            output: The parsed qsidd output.

        Returns:
            The joint and conditional probabilities.
        """
        competition_names = {
            transition.name: transition.names(competition=True)
            for transition in self.algorithm.transitions
        }
        certain = {
            name: (output.number(names.probability) or 0.0) >= 1.0
            for name, names in competition_names.items()
        }
        matrices: PairMatrices = {
            "joint_probability": {},
            "conditional_probability": {},
            "conditional_probability_given_absent": {},
        }

        for first in self.algorithm.transitions:
            for second in self.algorithm.transitions:
                if first is second:
                    continue

                pair = f"{first.letter} {second.letter}"
                given = f"{first.letter} given {second.letter}"
                given_not = f"{first.letter} given not {second.letter}"

                matrices["joint_probability"].setdefault(first.name, {})[
                    second.name
                ] = output.number(f"Joint probability {pair}")
                matrices["conditional_probability"].setdefault(first.name, {})[
                    second.name
                ] = output.number(f"Conditional probability {given}")
                matrices["conditional_probability_given_absent"].setdefault(
                    first.name, {}
                )[second.name] = (
                    None
                    if certain[second.name]
                    else output.number(f"Conditional probability {given_not}")
                )

        return matrices

    @staticmethod
    def _optional_number(output: QsiddOutput, key: str) -> float | None:
        """Read a value qsidd only prints for some algorithms.

        Args:
            output: The parsed qsidd output.
            key: Name of the value.

        Returns:
            The value, or None if qsidd did not print it.
        """
        return output.number(key) if key in output.values else None
