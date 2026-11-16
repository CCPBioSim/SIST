"""Run a SIST calculation with the qsidd binaries.

A run goes through these steps:

1. Look up the algorithm for the ``-a`` flag in :mod:`sist.transitions`.
2. For algorithms with cruciforms, score candidate inverted repeats with
   :mod:`sist.ir_finder`.
3. Run qsidd with the binary and flags the algorithm gives.
4. Write or print qsidd's text output, or, with ``--json``, build the
   structured result with :mod:`sist.results` and write it.

From Python, a run looks like::

    SistRunner(
        file="seq.fa",
        algorithm="M",
        temperature=310,
        superhelical_density=0.06,
        salt=0.01,
        energy_threshold=12,
        circular=False,
        nearest_neighbor=False,
        print_base_pair=False,
        print_parameters=False,
        print_ensemble_average=False,
        output_file=None,
        json_file="result.json",
    ).run()
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

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
        ir_string: The scored inverted repeats in the form ``qsidd -X``
            expects, such as ``"1,0,10000,|2982,25,17.27,|"``. Each entry is
            a start position, length and energy, ending in ``|``.
        irf_version: The IRF version, such as ``"3.09"``, if it could be read.
    """

    ir_string: str
    irf_version: str | None


class IRFinderScorer:
    """Score inverted repeats with the IR finder and IRF.

    IRF finds the candidate repeats in the sequence, and :class:`IRFinder`
    works out the energy of each one as a cruciform.

    Attributes:
        temperature: Temperature in kelvin.
        shape: ``"linear"`` or ``"circular"``.
    """

    def __init__(self, *, temperature: float, shape: str) -> None:
        """Set up the scorer with the conditions the energies are worked out at.

        Args:
            temperature: Temperature in kelvin.
            shape: ``"linear"`` or ``"circular"``.
        """
        self.temperature = temperature
        self.shape = shape

    def score(self, file: str) -> ScoredRepeats:
        """Score the inverted repeats in a sequence file.

        Args:
            file: Path to the FASTA sequence file.

        Returns:
            The scored inverted repeats, and the IRF version that found them.

        Raises:
            FileNotFoundError: If IRF is not installed.
        """
        finder = IRFinder(temperature=self.temperature, shape=self.shape)
        ir_string = finder.compute_cruciform_energy_string(file)

        return ScoredRepeats(ir_string, finder.irf_version)


class SubprocessQsidd:
    """Run qsidd as a subprocess.

    The exit code is returned rather than raised, so the runner can report a
    failed qsidd run as text or as a structured result.
    """

    def run(self, command: list[str]) -> subprocess.CompletedProcess[str]:
        """Run a qsidd command line and capture its output.

        Args:
            command: The command line to run, as built by
                :meth:`SistRunner._build_command`.

        Returns:
            The finished process, with its exit code and its stdout and stderr
            as text.
        """
        return subprocess.run(command, capture_output=True, text=True, check=False)


class SistRunner:
    """Run one SIST calculation and write or print its output.

    The attributes match the ``sist`` command's options. Text output works as
    it did before ``--json``: qsidd's output goes to ``output_file``, or to
    stdout. With ``json_file`` set, a structured result is written instead,
    and the text output only as well when ``output_file`` is also set.

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
        ir_scorer: IRFinderScorer | None = None,
        qsidd: SubprocessQsidd | None = None,
        result_writer: ResultWriter | None = None,
    ) -> None:
        """Set up the runner with one calculation's settings.

        The settings are described under the class attributes. The last three
        arguments set how inverted repeats are scored, how qsidd is run and
        how results are written. When they are None, an
        :class:`IRFinderScorer`, a :class:`SubprocessQsidd` and a
        :class:`ResultWriter` are used.
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
        """Molecule shape as the IR finder expects it.

        ``"circular"`` for a circular molecule, otherwise ``"linear"``.
        """
        return "circular" if self.circular else "linear"

    def run(self) -> None:
        """Run the calculation and write or print its output.

        Without a JSON file this writes qsidd's text output, exactly as before.
        With one, the structured result is written to it, and the text output
        is only written as well if ``output_file`` is set. The text output then
        comes from a second qsidd run with the flags the user chose, reusing
        the inverted repeats already scored.

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
        """Look up the algorithm for the ``-a`` flag in :data:`ALGORITHMS`.

        Returns:
            The algorithm, with its transitions, engine and qsidd flags.

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

        The output goes to ``output_file`` if it is set, otherwise to stdout.

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

        The steps are: read the sequence, score inverted repeats if needed,
        run qsidd with ``-b -p -r -P``, build the result and check it against
        the schema. A failure at any step, including a result that fails the
        check, becomes a failed result. That is written to ``json_file``
        before the error is raised again. Errors that are not a
        ``ResultError`` are recorded as ``INTERNAL_ERROR``.

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

        These fill in the ``calculation``, ``input`` and ``parameters``
        sections of the result.

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

        A missing IRF is recorded as ``IRF_NOT_FOUND`` in the failed result.

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
            ResultError: ``QSIDD_NOT_FOUND`` if the binary cannot be found, or
                ``QSIDD_FAILED`` with qsidd's exit code and stderr if it exits
                with an error.
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

        The parts come in this order: the binary, the boolean flags, the
        numeric parameters, the algorithm's own flags, the inverted repeats
        and the sequence file. For cruciform at full precision on a linear
        molecule this gives::

            trans_three -b -p -r -P -T 310.0 -s 0.06 -i 0.01 -t 12.0 -C \\
                -X <repeats> -f seq.fa

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

        For a structured result this is always ``-b -p -r -P``. For text output it
        is whichever of ``-b``, ``-p`` and ``-r`` were asked for. ``-n`` is added
        when nearest-neighbor energetics were asked for and the algorithm supports
        them, and ``-c`` for a circular molecule.

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
        """Build the qsidd numeric parameter flags shared by every algorithm.

        Returns:
            The temperature, superhelical density, salt concentration and
            energy threshold, as ``-T``, ``-s``, ``-i`` and ``-t`` with their
            values, such as ``["-T", "310.0", "-s", "0.06", "-i", "0.01", "-t",
            "12.0"]``.
        """
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
