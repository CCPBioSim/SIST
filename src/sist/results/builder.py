"""Put the parts of a result together for a successful or failed run."""

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
    """Build the result of one run.

    Create one builder when the run starts, so its start time is recorded,
    then call :meth:`read_sequence` and either :meth:`success` or
    :meth:`failure`.

    Attributes:
        settings: The run's settings.
        warnings: Warnings raised so far, in the order they were raised.
    """

    def __init__(self, settings: RunSettings, *, command: list[str] | None = None):
        """Start recording a run.

        Args:
            settings: The run's settings.
            command: The command line as given, or None when run from Python.
        """
        self.settings = settings
        self.warnings: list[ResultWarning] = []

        self._provenance = ProvenanceRecorder(
            engine=settings.algorithm.engine, command=command
        )
        self._reader = QsiddReader(settings.algorithm)

    def read_sequence(self) -> SequenceInput:
        """Read the input sequence and raise the warnings it and the settings call for.

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

        qsidd's warnings are recorded before checking for underflow, so a run
        that fails that check still keeps them.

        Args:
            output: The parsed qsidd output.
            sequence: The sequence that was analysed.
            ir_string: The scored inverted repeats passed to qsidd, if used.
            irf_version: The IRF version, if inverted repeats were scored.

        Returns:
            The result.

        Raises:
            ResultError: If the partition function underflowed, so qsidd's
                values are meaningless even though it exited normally.
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

        Everything known before the failure is kept, so a partition function
        underflow still records the derived parameters qsidd printed.

        Args:
            error: Why the run failed.
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
