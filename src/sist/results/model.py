"""The typed SIST result.

Each section of a result is a small dataclass with a ``to_dict`` method that
writes it in the shape of ``sist-result-v1.schema.json``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from sist.ir_finder import IRFinder
from sist.results.diagnostics import ResultError, ResultWarning
from sist.results.sequence import SequenceInput
from sist.transitions import Algorithm

SCHEMA_VERSION = "1.0.0"
"""Version of the result model this package writes."""

PairMatrices = dict[str, dict[str, dict[str, float | None]]]
"""Pairwise probabilities, read as ``matrices[name][first][second]``."""


@dataclass(frozen=True)
class RunSettings:
    """The settings for one SIST calculation.

    Attributes:
        file: Path to the input sequence file, as given.
        algorithm: The algorithm being run.
        temperature: Temperature in kelvin.
        superhelical_density: Superhelical density, as given.
        salt: Salt concentration in mol/L.
        energy_threshold: Energy threshold in kcal/mol.
        circular: Whether the molecule is circular.
        nearest_neighbor: Whether nearest-neighbor energetics were asked for.
    """

    file: str
    algorithm: Algorithm
    temperature: float
    superhelical_density: float
    salt: float
    energy_threshold: float
    circular: bool
    nearest_neighbor: bool

    @property
    def molecule_type(self) -> str:
        """``"circular"`` or ``"linear"``."""
        return "circular" if self.circular else "linear"

    def parameters(self) -> dict[str, Any]:
        """Write the settings as the calculation used them.

        The superhelical density is written as negative, as qsidd uses it.
        Runs that score inverted repeats also record the IR finder's fixed salt.

        Returns:
            The ``parameters`` section.
        """
        melting_energetics = None

        if "melting" in self.algorithm.transition_names:
            melting_energetics = (
                "nearest_neighbor" if self.nearest_neighbor else "copolymeric"
            )

        parameters: dict[str, Any] = {
            "temperature_k": self.temperature,
            "superhelical_density": -abs(self.superhelical_density),
            "salt_concentration_m": self.salt,
            "energy_threshold_kcal_per_mol": self.energy_threshold,
            "melting_energetics": melting_energetics,
        }

        if self.algorithm.needs_inverted_repeats:
            parameters["cruciform_salt_concentration_m"] = IRFinder.CRUCIFORM_SALT

        return parameters


@dataclass(frozen=True)
class Provenance:
    """How, where and when a result was made.

    Attributes:
        version: SIST version.
        git_sha: Git commit of the SIST source, if known.
        engine: qsidd binary used.
        irf_version: IRF version, for runs that score inverted repeats.
        python: Python version.
        platform: Operating system and architecture.
        command: The command line as given, or None when run from Python.
        started_at: When the run started.
        finished_at: When the run finished.
        wall_time_s: Elapsed time in seconds.
    """

    version: str
    git_sha: str | None
    engine: str
    irf_version: str | None
    python: str
    platform: str
    command: list[str] | None
    started_at: str
    finished_at: str
    wall_time_s: float

    def to_dict(self) -> dict[str, Any]:
        """Write the ``provenance`` section."""
        return {
            "software": {
                "name": "SIST",
                "version": self.version,
                "git_sha": self.git_sha,
                "engine": self.engine,
                "irf_version": self.irf_version,
            },
            "environment": {"python": self.python, "platform": self.platform},
            "command": self.command,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "wall_time_s": self.wall_time_s,
        }


@dataclass(frozen=True)
class InputDescription:
    """The input sequence, as analysed.

    Attributes:
        sequence_file: Path to the sequence file, as given.
        sequence_name: FASTA header, or None.
        sequence_length_bp: Base pairs analysed, or None if never read.
        analysed_sequence: The bases qsidd reports, or None if never read.
        sequence_sha256: SHA-256 of the analysed bases, or None if never read.
        molecule_type: ``"linear"`` or ``"circular"``.
    """

    sequence_file: str | None
    sequence_name: str | None
    sequence_length_bp: int | None
    analysed_sequence: str | None
    sequence_sha256: str | None
    molecule_type: str

    @classmethod
    def describe(
        cls, settings: RunSettings, sequence: SequenceInput | None
    ) -> InputDescription:
        """Describe the input, as far as it was read.

        Args:
            settings: The run's settings.
            sequence: The sequence, or None if it was never read.

        Returns:
            The input description.
        """
        if sequence is None:
            return cls(settings.file, None, None, None, None, settings.molecule_type)

        return cls(
            sequence_file=settings.file,
            sequence_name=sequence.name,
            sequence_length_bp=len(sequence.analysed),
            analysed_sequence=sequence.bases,
            sequence_sha256=sequence.sha256,
            molecule_type=settings.molecule_type,
        )

    def to_dict(self) -> dict[str, Any]:
        """Write the ``input`` section, leaving out a sequence that was never read."""
        section = asdict(self)

        if self.analysed_sequence is None:
            del section["analysed_sequence"]

        return section


@dataclass(frozen=True)
class DerivedParameters:
    """Model values qsidd works out from the settings.

    Attributes:
        linking_difference_turns: Imposed linking difference, in turns.
        torsional_stiffness_kcal_bp_per_rad2: Torsional stiffness.
        bat_kcal_per_mol: Copolymeric A-T melting energy, if printed.
        bgc_kcal_per_mol: Copolymeric G-C melting energy, if printed.
    """

    linking_difference_turns: float | None
    torsional_stiffness_kcal_bp_per_rad2: float | None
    bat_kcal_per_mol: float | None
    bgc_kcal_per_mol: float | None

    def to_dict(self) -> dict[str, Any]:
        """Write the ``derived_parameters`` section."""
        return asdict(self)


@dataclass(frozen=True)
class Summary:
    """Whole-molecule results.

    Attributes:
        ln_partition_function: Natural log of the partition function.
        partition_function: The partition function, or None if too small to
            store.
        average_number_of_runs: Expected total number of runs.
        states_by_number_of_runs: Entry k counts states with k + 1 runs.
        total_states: All states counted.
        pairwise: Pairwise probabilities, for competition runs.
    """

    ln_partition_function: float | None
    partition_function: float | None
    average_number_of_runs: float | None
    states_by_number_of_runs: list[int]
    total_states: int
    pairwise: PairMatrices | None = None

    def to_dict(self) -> dict[str, Any]:
        """Write the ``summary`` section."""
        summary: dict[str, Any] = {
            "ln_partition_function": self.ln_partition_function,
            "partition_function": self.partition_function,
            "average_number_of_runs": self.average_number_of_runs,
            "state_counts": {
                "by_number_of_runs": self.states_by_number_of_runs,
                "total": self.total_states,
            },
        }

        if self.pairwise is not None:
            summary["pairwise"] = self.pairwise

        return summary


@dataclass(frozen=True)
class TransitionResult:
    """Everything SIST reports about one transition.

    Attributes:
        nucleation_energy_kcal_per_mol: Nucleation energy.
        probability: Probability of at least one run.
        expected_runs: Expected number of runs.
        expected_transformed_bp: Expected number of base pairs in the transition.
        profile_probability: Probability at each position.
        free_energy_kcal_per_mol: G(x) at each position, where qsidd reports it.
    """

    nucleation_energy_kcal_per_mol: float | None
    probability: float | None
    expected_runs: float | None
    expected_transformed_bp: float | None
    profile_probability: list[float]
    free_energy_kcal_per_mol: list[float | None] | None = None

    def to_dict(self) -> dict[str, Any]:
        """Write one entry of the ``transitions`` section."""
        profile: dict[str, Any] = {"probability": self.profile_probability}

        if self.free_energy_kcal_per_mol is not None:
            profile["free_energy_kcal_per_mol"] = self.free_energy_kcal_per_mol

        return {
            "nucleation_energy_kcal_per_mol": self.nucleation_energy_kcal_per_mol,
            "probability": self.probability,
            "expected_runs": self.expected_runs,
            "expected_transformed_bp": self.expected_transformed_bp,
            "profile": profile,
        }


@dataclass(frozen=True)
class InvertedRepeat:
    """One candidate cruciform scored by the IR finder.

    Attributes:
        start: 1-based first position of the region.
        length_bp: Length of the region, arms and loop.
        energy_kcal_per_mol: Free energy to extrude the cruciform.
    """

    start: int
    length_bp: int
    energy_kcal_per_mol: float

    @classmethod
    def parse_all(cls, ir_string: str) -> list[InvertedRepeat]:
        """Read the scored inverted repeats from the string passed to qsidd.

        The string is ``start,length,energy,|`` repeated, and starts with a
        placeholder entry of length 0, which is skipped.

        Args:
            ir_string: The string the IR finder built for ``qsidd -X``.

        Returns:
            The scored inverted repeats.
        """
        inverted_repeats = []

        for entry in ir_string.split("|"):
            parts = entry.split(",")

            if len(parts) >= 3 and int(parts[1]) > 0:
                inverted_repeats.append(
                    cls(int(parts[0]), int(parts[1]), float(parts[2]))
                )

        return inverted_repeats

    def to_dict(self) -> dict[str, Any]:
        """Write one entry of the ``inverted_repeats`` section."""
        return asdict(self)


@dataclass(frozen=True)
class Diagnostics:
    """Values for checking a run.

    Attributes:
        scaling_factor: Offset qsidd adds to every Boltzmann exponent.
        min_energy_kcal_per_mol: Lowest state energy found.
        max_energy_kcal_per_mol: Highest state energy counted.
    """

    scaling_factor: float | None
    min_energy_kcal_per_mol: float | None
    max_energy_kcal_per_mol: float | None

    def to_dict(self) -> dict[str, Any]:
        """Write the ``diagnostics`` section."""
        return {
            "scaling_factor": self.scaling_factor,
            "energy_window": {
                "min_energy_kcal_per_mol": self.min_energy_kcal_per_mol,
                "max_energy_kcal_per_mol": self.max_energy_kcal_per_mol,
            },
        }


@dataclass(frozen=True)
class Result:
    """A complete SIST result, successful or failed.

    Attributes:
        provenance: How, where and when it was made.
        algorithm: The algorithm that was run.
        input: The input sequence.
        parameters: The settings as used, from :meth:`RunSettings.parameters`.
        warnings: Warnings raised during the run.
        derived_parameters: Model values, if known.
        summary: Whole-molecule results, for a successful run.
        transitions: Results for each transition, for a successful run.
        inverted_repeats: Scored inverted repeats, for runs that use them.
        diagnostics: Values for checking a successful run.
        error: Why the run failed, for a failed run.
    """

    provenance: Provenance
    algorithm: Algorithm
    input: InputDescription
    parameters: dict[str, Any]
    warnings: list[ResultWarning]
    derived_parameters: DerivedParameters | None = None
    summary: Summary | None = None
    transitions: dict[str, TransitionResult] | None = None
    inverted_repeats: list[InvertedRepeat] | None = None
    diagnostics: Diagnostics | None = None
    error: ResultError | None = None

    @property
    def status(self) -> str:
        """``"failed"`` if there is an error, otherwise ``"success"``."""
        return "failed" if self.error is not None else "success"

    def to_dict(self) -> dict[str, Any]:
        """Write the result in the shape of ``sist-result-v1.schema.json``.

        Returns:
            The result, with its sections in the documented order.
        """
        result: dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            "status": self.status,
            "provenance": self.provenance.to_dict(),
            "calculation": {
                "algorithm": self.algorithm.flag,
                "transitions": list(self.algorithm.transition_names),
            },
            "input": self.input.to_dict(),
            "parameters": self.parameters,
            "derived_parameters": _to_dict_or_none(self.derived_parameters),
            "summary": _to_dict_or_none(self.summary),
            "transitions": (
                {name: result.to_dict() for name, result in self.transitions.items()}
                if self.transitions is not None
                else None
            ),
        }

        if self.inverted_repeats is not None:
            result["inverted_repeats"] = [
                inverted_repeat.to_dict() for inverted_repeat in self.inverted_repeats
            ]

        if self.diagnostics is not None:
            result["diagnostics"] = self.diagnostics.to_dict()

        result["warnings"] = [warning.to_dict() for warning in self.warnings]
        result["error"] = _to_dict_or_none(self.error)

        return result


def _to_dict_or_none(section: Any) -> dict[str, Any] | None:
    """Write a section with its own ``to_dict``, or None if it is missing.

    Args:
        section: A section object with a ``to_dict`` method, or None.

    Returns:
        The section as a dictionary, or None.
    """
    return section.to_dict() if section is not None else None
