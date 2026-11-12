"""Read the input sequence the way qsidd does."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from sist.results.diagnostics import ErrorCode, ResultError


@dataclass(frozen=True)
class SequenceInput:
    """The sequence read from the input file, as qsidd will analyse it.

    Header lines are skipped, A, C, G, T and N are kept in any case, and every
    other character is dropped, as qsidd does.

    Attributes:
        name: The FASTA header without the leading ``>``, or None if the file
            has no header.
        analysed: The bases qsidd analyses, in upper case.
        dropped: Number of characters dropped that are not A, C, G, T or N.
    """

    name: str | None
    analysed: str
    dropped: int

    @classmethod
    def read(cls, path: Path) -> SequenceInput:
        """Read a sequence file the way qsidd does.

        Args:
            path: Path to the sequence file.

        Returns:
            The sequence that qsidd will analyse.

        Raises:
            ResultError: If the file does not exist, holds more than one FASTA
                record, or has no bases to analyse.
        """
        if not path.is_file():
            raise ResultError(
                ErrorCode.INPUT_NOT_FOUND, f"The sequence file {path} does not exist."
            )

        lines = path.read_text(encoding="utf-8").splitlines()
        headers = [line for line in lines if line.startswith(">")]

        if len(headers) > 1:
            raise ResultError(
                ErrorCode.MULTIPLE_SEQUENCES,
                f"The sequence file has {len(headers)} FASTA records. SIST analyses "
                "one sequence per run.",
                details={"records": len(headers)},
            )

        letters = "".join(
            "".join(line.split()) for line in lines if not line.startswith(">")
        ).upper()
        analysed = "".join(base for base in letters if base in "ACGTN")

        if not analysed:
            raise ResultError(
                ErrorCode.EMPTY_SEQUENCE,
                f"The sequence file {path} has no bases to analyse.",
            )

        name = headers[0][1:].strip() if headers else None

        return cls(name=name, analysed=analysed, dropped=len(letters) - len(analysed))

    @property
    def bases(self) -> str:
        """The analysed sequence with N replaced by G, as qsidd reports it."""
        return self.analysed.replace("N", "G")

    @property
    def sha256(self) -> str:
        """SHA-256 of the bases qsidd reports, for comparing results."""
        return hashlib.sha256(self.bases.encode("ascii")).hexdigest()

    @property
    def n_positions(self) -> list[int]:
        """1-based positions of every N in the analysed sequence."""
        return [index + 1 for index, base in enumerate(self.analysed) if base == "N"]
