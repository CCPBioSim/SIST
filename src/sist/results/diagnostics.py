"""Warning and error codes, and the checks that raise warnings.

Every code a result can hold is listed here once, along with the part of SIST
that raises it, so codes cannot drift between the code, the tests and the docs.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import TYPE_CHECKING, Any

from sist.ir_finder import IRFinder

if TYPE_CHECKING:
    from sist.results.sequence import SequenceInput
    from sist.transitions import Algorithm


class WarningCode(StrEnum):
    """Every warning a result can hold."""

    THRESHOLD_TOO_LOW = "threshold_too_low"
    THRESHOLD_TOO_HIGH = "threshold_too_high"
    SUPERHELICAL_DENSITY_OUT_OF_RANGE = "superhelical_density_out_of_range"
    TEMPERATURE_OUT_OF_RANGE = "temperature_out_of_range"
    SALT_CONCENTRATION_OUT_OF_RANGE = "salt_concentration_out_of_range"
    SEQUENCE_TOO_SHORT = "sequence_too_short"
    SEQUENCE_TOO_LONG = "sequence_too_long"
    AMBIGUOUS_BASES_REPLACED = "ambiguous_bases_replaced"
    UNRECOGNISED_CHARACTERS_IGNORED = "unrecognised_characters_ignored"
    NEAREST_NEIGHBOR_IGNORED = "nearest_neighbor_ignored"
    CRUCIFORM_SALT_FIXED = "cruciform_salt_fixed"
    NO_INVERTED_REPEATS_FOUND = "no_inverted_repeats_found"

    @property
    def source(self) -> str:
        """Part of SIST that raises this warning."""
        return WARNING_SOURCES.get(self, "sist")


class ErrorCode(StrEnum):
    """Every reason a run can fail."""

    INPUT_NOT_FOUND = "input_not_found"
    EMPTY_SEQUENCE = "empty_sequence"
    MULTIPLE_SEQUENCES = "multiple_sequences"
    IRF_NOT_FOUND = "irf_not_found"
    QSIDD_NOT_FOUND = "qsidd_not_found"
    QSIDD_FAILED = "qsidd_failed"
    PARTITION_FUNCTION_UNDERFLOW = "partition_function_underflow"
    INTERNAL_ERROR = "internal_error"

    @property
    def stage(self) -> str:
        """Step of the run that fails with this code."""
        return ERROR_STAGES[self]


WARNING_SOURCES: dict[WarningCode, str] = {
    WarningCode.THRESHOLD_TOO_LOW: "qsidd",
    WarningCode.THRESHOLD_TOO_HIGH: "qsidd",
    WarningCode.SUPERHELICAL_DENSITY_OUT_OF_RANGE: "qsidd",
    WarningCode.TEMPERATURE_OUT_OF_RANGE: "qsidd",
    WarningCode.SALT_CONCENTRATION_OUT_OF_RANGE: "qsidd",
    WarningCode.NO_INVERTED_REPEATS_FOUND: "ir_finder",
}

ERROR_STAGES: dict[ErrorCode, str] = {
    ErrorCode.INPUT_NOT_FOUND: "input",
    ErrorCode.EMPTY_SEQUENCE: "input",
    ErrorCode.MULTIPLE_SEQUENCES: "input",
    ErrorCode.IRF_NOT_FOUND: "ir_finder",
    ErrorCode.QSIDD_NOT_FOUND: "qsidd",
    ErrorCode.QSIDD_FAILED: "qsidd",
    ErrorCode.PARTITION_FUNCTION_UNDERFLOW: "qsidd",
    ErrorCode.INTERNAL_ERROR: "internal",
}

QSIDD_WARNINGS: dict[str, WarningCode] = {
    "threshold is too small, results may be inaccurate": WarningCode.THRESHOLD_TOO_LOW,
    "threshold is too high, execution time may be very long": (
        WarningCode.THRESHOLD_TOO_HIGH
    ),
    "superhelical density is outside of physiological range": (
        WarningCode.SUPERHELICAL_DENSITY_OUT_OF_RANGE
    ),
    "temperature is outside of physiological range": (
        WarningCode.TEMPERATURE_OUT_OF_RANGE
    ),
    "salt concentration is outside of physiological range": (
        WarningCode.SALT_CONCENTRATION_OUT_OF_RANGE
    ),
}
"""qsidd's warning text and its code. qsidd's two sequence length warnings are
left out, because it checks the size of the file rather than the number of base
pairs, so :func:`check_sequence` raises those instead."""

MIN_SEQUENCE_LENGTH = 1500
"""Shortest recommended sequence, in base pairs."""

MAX_SEQUENCE_LENGTH = 10000
"""Longest recommended sequence, in base pairs."""


@dataclass(frozen=True)
class ResultWarning:
    """A warning recorded in a result.

    Attributes:
        code: What kind of warning it is.
        message: What it means, written for people.
        details: Extra values about the warning.
    """

    code: WarningCode
    message: str
    details: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        """Write the warning in the shape the result schema expects.

        Returns:
            One entry of a result's ``warnings`` list.
        """
        warning: dict[str, Any] = {
            "code": str(self.code),
            "message": self.message,
            "source": self.code.source,
        }

        if self.details:
            warning["details"] = self.details

        return warning


class ResultError(Exception):
    """A run that could not produce a result, and why.

    Attributes:
        code: Why the run failed.
        message: What went wrong, written for people.
        exit_code: Exit code of the failing external program, if any.
        stderr: Standard error of the failing external program, if any.
        details: Extra values about the error.
    """

    def __init__(
        self,
        code: ErrorCode,
        message: str,
        *,
        exit_code: int | None = None,
        stderr: str | None = None,
        details: dict[str, Any] | None = None,
    ) -> None:
        """Record why the run failed.

        Args:
            code: Why the run failed.
            message: What went wrong.
            exit_code: Exit code of the failing external program, if any.
            stderr: Standard error of the failing external program, if any.
            details: Extra values about the error.
        """
        super().__init__(message)

        self.code = code
        self.message = message
        self.exit_code = exit_code
        self.stderr = stderr
        self.details = details

    def to_dict(self) -> dict[str, Any]:
        """Write the error in the shape the result schema expects.

        Returns:
            The ``error`` section of a failed result.
        """
        error: dict[str, Any] = {
            "code": str(self.code),
            "message": self.message,
            "stage": self.code.stage,
            "exit_code": self.exit_code,
            "stderr": self.stderr,
        }

        if self.details:
            error["details"] = self.details

        return error


def check_sequence(sequence: SequenceInput) -> list[ResultWarning]:
    """Raise the warnings the input sequence calls for.

    Args:
        sequence: The sequence qsidd will analyse.

    Returns:
        Warnings about ambiguous bases, dropped characters and the length.
    """
    warnings = []
    length = len(sequence.analysed)
    n_positions = sequence.n_positions

    if n_positions:
        warnings.append(
            ResultWarning(
                WarningCode.AMBIGUOUS_BASES_REPLACED,
                f"{len(n_positions)} N base(s) in the sequence were analysed as G.",
                {"count": len(n_positions), "positions": n_positions},
            )
        )

    if sequence.dropped:
        warnings.append(
            ResultWarning(
                WarningCode.UNRECOGNISED_CHARACTERS_IGNORED,
                f"{sequence.dropped} character(s) other than A, C, G, T and N were "
                "dropped, so later positions are shifted.",
                {"count": sequence.dropped},
            )
        )

    if length < MIN_SEQUENCE_LENGTH:
        warnings.append(
            ResultWarning(
                WarningCode.SEQUENCE_TOO_SHORT,
                f"The sequence is {length} bp, shorter than the recommended "
                f"{MIN_SEQUENCE_LENGTH} bp.",
                {"sequence_length_bp": length},
            )
        )

    if length > MAX_SEQUENCE_LENGTH:
        warnings.append(
            ResultWarning(
                WarningCode.SEQUENCE_TOO_LONG,
                f"The sequence is {length} bp, longer than the recommended "
                f"{MAX_SEQUENCE_LENGTH} bp.",
                {"sequence_length_bp": length},
            )
        )

    return warnings


def check_settings(
    algorithm: Algorithm, *, nearest_neighbor: bool, salt: float
) -> list[ResultWarning]:
    """Raise the warnings for settings the algorithm cannot use as given.

    Args:
        algorithm: The algorithm being run.
        nearest_neighbor: Whether nearest-neighbor energetics were asked for.
        salt: Salt concentration in mol/L.

    Returns:
        Warnings about ignored energetics and the fixed cruciform salt.
    """
    warnings = []

    if nearest_neighbor and not algorithm.supports_nearest_neighbor:
        warnings.append(
            ResultWarning(
                WarningCode.NEAREST_NEIGHBOR_IGNORED,
                "Nearest-neighbor energetics only apply to melting, so they were "
                "ignored.",
            )
        )

    if algorithm.needs_inverted_repeats and salt != IRFinder.CRUCIFORM_SALT:
        warnings.append(
            ResultWarning(
                WarningCode.CRUCIFORM_SALT_FIXED,
                f"Inverted repeats are always scored at {IRFinder.CRUCIFORM_SALT} M "
                f"salt, not the {salt} M used by qsidd.",
            )
        )

    return warnings


def qsidd_warnings(messages: list[str]) -> list[ResultWarning]:
    """Turn qsidd's warning lines into coded warnings.

    Args:
        messages: The text of each ``WARNING:`` line qsidd printed.

    Returns:
        A warning for each line with a known code, in the same order.
    """
    return [
        ResultWarning(QSIDD_WARNINGS[message], f"{message[0].upper()}{message[1:]}.")
        for message in messages
        if message in QSIDD_WARNINGS
    ]
