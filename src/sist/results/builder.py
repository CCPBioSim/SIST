"""Put the parts of a result together for a successful or failed run.

:class:`ResultBuilder` brings together the other parts of :mod:`sist.results`:

- :class:`SequenceInput` reads the input sequence.
- :mod:`sist.results.diagnostics` checks the sequence, the settings and
  qsidd's output, and raises warnings.
- :class:`QsiddReader` reads the values out of the parsed qsidd output.
- :class:`ProvenanceRecorder` records how, where and when the run was made.

The result it returns is a :class:`Result`, ready for :class:`ResultWriter`.
"""

from __future__ import annotations

from pathlib import Path

from sist.results.diagnostics import (
    ErrorCode,
    ResultError,
    ResultWarning,
    WarningCode,
    check_sequence,
    check_settings,
    qsidd_warnings,
)
from sist.results.model import (
    DerivedParameters,
    InputDescription,
    InvertedRepeat,
    Result,
    RunSettings,
)
from sist.results.provenance import ProvenanceRecorder
from sist.results.qsidd_output import QsiddOutput
from sist.results.reader import QsiddReader
from sist.results.sequence import SequenceInput


class ResultBuilder:
    """Build the result of one run, successful or failed.

    A run uses one builder from start to finish::

        builder = ResultBuilder(settings, command=argv)
        sequence = builder.read_sequence()
        result = builder.success(output, sequence, ir_string=ir_string)

    ``read_sequence`` and ``success`` add warnings to ``warnings`` as they go.
    If any step raises a ``ResultError``, ``failure`` builds the failed result
    instead, keeping the warnings and anything already read.

    Attributes:
        settings: The run's settings.
        warnings: Warnings raised so far, in the order they were raised.
    """

    def __init__(self, settings: RunSettings, *, command: list[str] | None = None):
        """Start recording a run, and its start time.

        Args:
            settings: The run's settings, which also pick the qsidd binary
                and how its output is read.
            command: The command line as given, such as
                ``["sist", "-a", "M", "-f", "seq.fa", "--json", "out.json"]``,
                or None when run from Python.
        """
        self.settings = settings
        self.warnings: list[ResultWarning] = []

        self._provenance = ProvenanceRecorder(
            engine=settings.algorithm.engine, command=command
        )
        self._reader = QsiddReader(settings.algorithm)

    def read_sequence(self) -> SequenceInput:
        """Read the input sequence and raise the warnings it and the settings call for.

        The sequence warnings come first, such as characters dropped or ``N``
        read as ``G``, then the settings warnings, such as ``-n`` given for an
        algorithm that ignores it.

        Returns:
            The sequence that qsidd will analyse.

        Raises:
            ResultError: If the sequence cannot be analysed.
        """
        sequence = SequenceInput.read(Path(self.settings.file))

        self.warnings.extend(check_sequence(sequence))
        self.warnings.extend(
            check_settings(
                self.settings.algorithm,
                nearest_neighbor=self.settings.nearest_neighbor,
                salt=self.settings.salt,
            )
        )

        return sequence

    def success(
        self,
        output: QsiddOutput,
        sequence: SequenceInput,
        *,
        ir_string: str | None = None,
        irf_version: str | None = None,
    ) -> Result:
        """Build the result of a run that finished.

        Adds qsidd's own warnings, such as an energy threshold that is too
        low, then checks the partition function did not underflow. For
        algorithms with cruciforms, the inverted repeats are parsed from
        ``ir_string``, and a warning is added if there are none.

        Args:
            output: The parsed qsidd output, from a run with ``-b -p -r -P``.
            sequence: The sequence that was analysed.
            ir_string: The scored inverted repeats passed to qsidd with
                ``-X``, if used.
            irf_version: The IRF version, such as ``"3.09"``, if inverted
                repeats were scored.

        Returns:
            The successful result, with every section filled in.

        Raises:
            ResultError: If the partition function underflowed to zero.
        """
        self.warnings.extend(qsidd_warnings(output.warnings))

        if self._reader.underflowed(output):
            raise ResultError(
                ErrorCode.PARTITION_FUNCTION_UNDERFLOW,
                "The partition function underflowed to zero, so qsidd could not "
                "work out any probabilities. Try a lower superhelical density or "
                "a shorter sequence.",
                exit_code=0,
                details={"min_energy_kcal_per_mol": output.number("MinE")},
            )

        inverted_repeats = None

        if self.settings.algorithm.needs_inverted_repeats:
            inverted_repeats = InvertedRepeat.parse_all(ir_string or "")

            if not inverted_repeats:
                self.warnings.append(
                    ResultWarning(
                        WarningCode.NO_INVERTED_REPEATS_FOUND,
                        "IRF found no candidate inverted repeats, so no cruciforms "
                        "can form.",
                    )
                )

        return Result(
            provenance=self._provenance.finish(irf_version=irf_version),
            algorithm=self.settings.algorithm,
            input=InputDescription.describe(self.settings, sequence),
            parameters=self.settings.parameters(),
            warnings=self.warnings,
            derived_parameters=self._reader.derived_parameters(output),
            summary=self._reader.summary(output),
            transitions=self._reader.transitions(output),
            inverted_repeats=inverted_repeats,
            diagnostics=self._reader.diagnostics(output),
        )

    def failure(
        self,
        error: ResultError,
        *,
        sequence: SequenceInput | None = None,
        output: QsiddOutput | None = None,
        irf_version: str | None = None,
    ) -> Result:
        """Build the result of a run that failed.

        Keeps everything known before the failure: the settings, the
        warnings so far, the input if it was read, and any derived parameters
        qsidd printed. The summary, transitions and diagnostics are left out,
        and the error is recorded in their place.

        Args:
            error: Why the run failed, with its code, stage and message.
            sequence: The sequence, if it was read before the failure.
            output: The parsed qsidd output, if qsidd ran.
            irf_version: The IRF version, if IRF ran.

        Returns:
            The failed result.
        """
        return Result(
            provenance=self._provenance.finish(irf_version=irf_version),
            algorithm=self.settings.algorithm,
            input=InputDescription.describe(self.settings, sequence),
            parameters=self.settings.parameters(),
            warnings=self.warnings,
            derived_parameters=self._known_derived_parameters(output),
            error=error,
        )

    def _known_derived_parameters(
        self, output: QsiddOutput | None
    ) -> DerivedParameters | None:
        """Read the derived parameters for a failed run, if they can be read.

        qsidd may have stopped before printing them, or printed only some, in
        which case None is returned.

        Args:
            output: The parsed qsidd output, or None if qsidd did not run.

        Returns:
            The derived parameters, or None if they are not known.
        """
        if output is None:
            return None

        try:
            return self._reader.derived_parameters(output)
        except (KeyError, ValueError):
            return None
