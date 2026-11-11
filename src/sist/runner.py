"""Orchestration that dispatches a SIST calculation to the qsidd binaries.

Algorithms that need inverted repeats first score candidates with
:mod:`sist.ir_finder`. qsidd is then run with the binary and flags
:mod:`sist.transitions` gives for the algorithm.

The text output is written as before. When a JSON file is asked for, qsidd is
run at full precision and the result is built by :mod:`sist.results`.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from sist.binaries import BinaryLocator
from sist.ir_finder import IRFinder
from sist.results.builder import ResultBuilder
from sist.results.diagnostics import ErrorCode, ResultError
from sist.results.model import RunSettings
from sist.results.qsidd_output import QsiddOutput
from sist.results.sequence import SequenceInput
from sist.results.writer import ResultWriter
from sist.transitions import ALGORITHMS, Algorithm


@dataclass(frozen=True)
class ScoredRepeats:
    """Inverted repeats scored for qsidd, and the IRF version that found them.

    Attributes:
        ir_string: The scored inverted repeats in the form ``qsidd -X`` expects.
        irf_version: The IRF version, if it could be read.
    """

    ir_string: str
    irf_version: str | None


class InvertedRepeatScorer(Protocol):
    """Something that scores candidate inverted repeats in a sequence file."""

    def score(self, file: str) -> ScoredRepeats:
        """Score the inverted repeats in a sequence file.

        Args:
            file: Path to the sequence file.

        Returns:
            The scored inverted repeats.
        """
        ...


class QsiddExecutor(Protocol):
    """Something that runs a qsidd command line."""

    def run(self, command: list[str]) -> subprocess.CompletedProcess[str]:
        """Run a qsidd command line.

        Args:
            command: The command line to run.

        Returns:
            The finished process, with its output as text.
        """
        ...


class IRFinderScorer:
    """Score inverted repeats with the IR finder and IRF.

    Attributes:
        temperature: Temperature in kelvin.
        shape: ``"linear"`` or ``"circular"``.
    """

    def __init__(self, *, temperature: float, shape: str) -> None:
        """Set up the scorer.

        Args:
            temperature: Temperature in kelvin.
            shape: ``"linear"`` or ``"circular"``.
        """
        self.temperature = temperature
        self.shape = shape

    def score(self, file: str) -> ScoredRepeats:
        """Score the inverted repeats in a sequence file.

        Args:
            file: Path to the sequence file.

        Returns:
            The scored inverted repeats.
        """
        finder = IRFinder(temperature=self.temperature, shape=self.shape)
        ir_string = finder.compute_cruciform_energy_string(file)

        return ScoredRepeats(ir_string, finder.irf_version)


class SubprocessQsidd:
    """Run qsidd as a subprocess."""

    def run(self, command: list[str]) -> subprocess.CompletedProcess[str]:
        """Run a qsidd command line and capture its output.

        Args:
            command: The command line to run.

        Returns:
            The finished process, with its output as text.
        """
        return subprocess.run(command, capture_output=True, text=True, check=False)


class SistRunner:
    """Coordinate a single SIST calculation: build, run, and report on it.

    Attributes:
        file: Path to the input FASTA sequence file.
        algorithm: One of ``"M"`` (melting), ``"Z"`` (Z-DNA),
            ``"C"`` (cruciform), or ``"A"`` (competition).
        temperature: Temperature in Kelvin.
        superhelical_density: Superhelical density (sigma).
        salt: Salt (ionic strength) concentration in mol/L.
        energy_threshold: Energy threshold for the transition search.
        circular: Whether the input molecule is circular (default linear).
        nearest_neighbor: Use nearest-neighbor melting energetics.
            Ignored for the Z-DNA and cruciform algorithms.
        print_base_pair: Print the base pair for each position.
        print_parameters: Print algorithm parameters.
        print_ensemble_average: Print ensemble average results.
        output_file: Path to write output to, or None to print to stdout.
        json_file: Path to write the structured result to, or None to skip it.
        command: The command line as given, recorded in the structured result.
    """

    def __init__(
        self,
        *,
        file: str,
        algorithm: str,
        temperature: float,
        superhelical_density: float,
        salt: float,
        energy_threshold: float,
        circular: bool,
        nearest_neighbor: bool,
        print_base_pair: bool,
        print_parameters: bool,
        print_ensemble_average: bool,
        output_file: str | None,
        json_file: str | None = None,
        command: list[str] | None = None,
        ir_scorer: InvertedRepeatScorer | None = None,
        qsidd: QsiddExecutor | None = None,
        result_writer: ResultWriter | None = None,
    ) -> None:
        """Initialise the runner with a single calculation's configuration.

        The last three arguments replace how inverted repeats are scored, how
        qsidd is run and how results are written.
        """
        self.file = file
        self.algorithm = algorithm
        self.temperature = temperature
        self.superhelical_density = superhelical_density
        self.salt = salt
        self.energy_threshold = energy_threshold
        self.circular = circular
        self.nearest_neighbor = nearest_neighbor
        self.print_base_pair = print_base_pair
        self.print_parameters = print_parameters
        self.print_ensemble_average = print_ensemble_average
        self.output_file = output_file
        self.json_file = json_file
        self.command = command

        self._binaries = BinaryLocator()
        self._ir_scorer = ir_scorer or IRFinderScorer(
            temperature=temperature, shape=self.shape
        )
        self._qsidd = qsidd or SubprocessQsidd()
        self._result_writer = result_writer or ResultWriter()

    @property
    def shape(self) -> str:
        """Molecule shape as the qsidd binaries expect it: linear or circular."""
        return "circular" if self.circular else "linear"

    def run(self) -> None:
        """Run the calculation and write or print its output.

        Without a JSON file this writes qsidd's text output, exactly as before.
        With one, the structured result is written to it, and the text output
        is only written as well if ``output_file`` is set. The text output then
        comes from a second qsidd run without full precision.

        Raises:
            ValueError: If ``algorithm`` is not one of M, Z, C, or A.
            subprocess.CalledProcessError: If the qsidd binary exits
                non-zero while writing text output.
            ResultError: If the run failed while writing a structured result.
                The failed result is written before this is raised.
        """
        algorithm = self._algorithm()

        if self.json_file:
            repeats = self._run_json(algorithm)

            if self.output_file:
                self._write_text(algorithm, repeats)

            return

        repeats = (
            self._ir_scorer.score(self.file)
            if algorithm.needs_inverted_repeats
            else None
        )
        self._write_text(algorithm, repeats)

    def _algorithm(self) -> Algorithm:
        """Look up the algorithm to run.

        Returns:
            The algorithm.

        Raises:
            ValueError: If the algorithm flag is not known.
        """
        algorithm = ALGORITHMS.get(self.algorithm)

        if algorithm is None:
            raise ValueError(
                f"Unknown algorithm: {self.algorithm!r}. Expected one of "
                f"{', '.join(ALGORITHMS)}."
            )

        return algorithm

    def _write_text(self, algorithm: Algorithm, repeats: ScoredRepeats | None) -> None:
        """Run qsidd as usual and write or print its text output.

        Args:
            algorithm: The algorithm to run.
            repeats: The scored inverted repeats, if the algorithm uses them.

        Raises:
            subprocess.CalledProcessError: If qsidd exits non-zero.
        """
        command = self._build_command(algorithm, repeats, precise=False)
        completed = self._qsidd.run(command)

        if completed.returncode != 0:
            raise subprocess.CalledProcessError(
                completed.returncode, command, completed.stdout, completed.stderr
            )

        if self.output_file:
            Path(self.output_file).write_text(completed.stdout, encoding="utf-8")
        else:
            print(completed.stdout, end="")

    def _run_json(self, algorithm: Algorithm) -> ScoredRepeats | None:
        """Run qsidd at full precision and write the structured result.

        Turns every failure, including a result that fails validation, into a
        failed result and writes it before raising the error again.

        Args:
            algorithm: The algorithm to run.

        Returns:
            The scored inverted repeats, for reuse by a text run.

        Raises:
            ResultError: If the run failed.
        """
        builder = ResultBuilder(self._settings(algorithm), command=self.command)
        sequence: SequenceInput | None = None
        output: QsiddOutput | None = None
        repeats: ScoredRepeats | None = None

        try:
            sequence = builder.read_sequence()
            repeats = self._score_for_json(algorithm)
            output = self._run_precise_qsidd(algorithm, repeats)
            result = builder.success(
                output,
                sequence,
                ir_string=repeats.ir_string if repeats else None,
                irf_version=repeats.irf_version if repeats else None,
            )
            self._result_writer.check(result)
        except Exception as error:
            failure = (
                error
                if isinstance(error, ResultError)
                else ResultError(ErrorCode.INTERNAL_ERROR, str(error))
            )
            self._result_writer.write(
                builder.failure(
                    failure,
                    sequence=sequence,
                    output=output,
                    irf_version=repeats.irf_version if repeats else None,
                ),
                Path(str(self.json_file)),
            )

            if failure is error:
                raise

            raise failure from error

        self._result_writer.write(result, Path(str(self.json_file)))

        return repeats

    def _settings(self, algorithm: Algorithm) -> RunSettings:
        """Collect the settings a structured result records.

        Args:
            algorithm: The algorithm being run.

        Returns:
            The run's settings.
        """
        return RunSettings(
            file=self.file,
            algorithm=algorithm,
            temperature=self.temperature,
            superhelical_density=self.superhelical_density,
            salt=self.salt,
            energy_threshold=self.energy_threshold,
            circular=self.circular,
            nearest_neighbor=self.nearest_neighbor,
        )

    def _score_for_json(self, algorithm: Algorithm) -> ScoredRepeats | None:
        """Score inverted repeats, reporting a missing IRF as a result error.

        Args:
            algorithm: The algorithm being run.

        Returns:
            The scored inverted repeats, or None if the algorithm does not use
            them.

        Raises:
            ResultError: If IRF is not installed.
        """
        if not algorithm.needs_inverted_repeats:
            return None

        try:
            return self._ir_scorer.score(self.file)
        except FileNotFoundError as error:
            raise ResultError(
                ErrorCode.IRF_NOT_FOUND,
                f"The Inverted Repeats Finder executable '{IRFinder.IRF_EXECUTABLE}' "
                "was not found on PATH. Algorithms that need inverted repeats "
                "cannot run without it.",
            ) from error

    def _run_precise_qsidd(
        self, algorithm: Algorithm, repeats: ScoredRepeats | None
    ) -> QsiddOutput:
        """Run qsidd with every output option and full precision.

        Args:
            algorithm: The algorithm to run.
            repeats: The scored inverted repeats, if the algorithm uses them.

        Returns:
            The parsed qsidd output.

        Raises:
            ResultError: If qsidd cannot be found or exits with an error.
        """
        try:
            command = self._build_command(algorithm, repeats, precise=True)
        except FileNotFoundError as error:
            raise ResultError(ErrorCode.QSIDD_NOT_FOUND, str(error)) from error

        completed = self._qsidd.run(command)

        if completed.returncode != 0:
            raise ResultError(
                ErrorCode.QSIDD_FAILED,
                f"qsidd exited with code {completed.returncode}.",
                exit_code=completed.returncode,
                stderr=completed.stderr or None,
            )

        return QsiddOutput.parse(completed.stdout)

    def _build_command(
        self, algorithm: Algorithm, repeats: ScoredRepeats | None, *, precise: bool
    ) -> list[str]:
        """Build the qsidd command line for an algorithm.

        Args:
            algorithm: The algorithm to run.
            repeats: The scored inverted repeats, if the algorithm uses them.
            precise: Ask for every output section at full precision, in place
                of the flags the user chose.

        Returns:
            The qsidd command line.

        Raises:
            FileNotFoundError: If the qsidd binary cannot be found.
        """
        command = [
            str(self._binaries.binary(algorithm.engine)),
            *self._flags(algorithm, precise=precise),
            *self._shared_parameters(),
            *algorithm.qsidd_flags,
        ]

        if algorithm.needs_inverted_repeats:
            command.extend(["-X", repeats.ir_string if repeats else ""])

        return [*command, "-f", self.file]

    def _flags(self, algorithm: Algorithm, *, precise: bool) -> list[str]:
        """Build the qsidd boolean flags.

        Args:
            algorithm: The algorithm being run.
            precise: Ask for every output section at full precision.

        Returns:
            The flags.
        """
        flags = []

        if precise:
            flags.extend(["-b", "-p", "-r", "-P"])
        else:
            if self.print_base_pair:
                flags.append("-b")
            if self.print_parameters:
                flags.append("-p")
            if self.print_ensemble_average:
                flags.append("-r")

        if self.nearest_neighbor and algorithm.supports_nearest_neighbor:
            flags.append("-n")
        if self.circular:
            flags.append("-c")

        return flags

    def _shared_parameters(self) -> list[str]:
        """Build the qsidd numeric parameter flags shared by every algorithm."""
        return [
            "-T",
            str(self.temperature),
            "-s",
            str(self.superhelical_density),
            "-i",
            str(self.salt),
            "-t",
            str(self.energy_threshold),
        ]
