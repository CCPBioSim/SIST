"""Read the input sequence the way qsidd does.

:meth:`SequenceInput.read` turns a FASTA file into the bases qsidd will analyse,
and raises a :class:`ResultError` for input that cannot be analysed.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from sist.results.diagnostics import ErrorCode, ResultError


@dataclass(frozen=True)
class SequenceInput:
    """The sequence read from the input file, as qsidd will analyse it.

    Header lines are skipped, A, C, G, T and N are kept in any case, and every
    other character is dropped, as qsidd does. For a file holding
    ``">example"`` then ``"acgRNt"``:

    >>> sequence = SequenceInput(name="example", analysed="ACGNT", dropped=1)
    >>> sequence.bases
    'ACGGT'
    >>> sequence.n_positions
    [4]

    Attributes:
        name: The FASTA header without the leading ``>``, or None if the file
            has no header.
        analysed: The bases qsidd analyses, in upper case, with any ``N`` kept.
        dropped: Number of characters dropped that are not A, C, G, T or N.
    """

    name: str | None
    analysed: str
    dropped: int

    @classmethod
    def read(cls, path: Path) -> SequenceInput:
        """Read a sequence file the way qsidd does.

        Line breaks and other whitespace are ignored, so a sequence split over
        many lines is read as one.

        Args:
            path: Path to the sequence file.

        Returns:
            The sequence that qsidd will analyse.

        Raises:
            ResultError: With ``input_not_found`` if the file does not exist,
                ``multiple_sequences`` if it holds more than one FASTA record,
                or ``empty_sequence`` if it has no bases to analyse.
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
        """The analysed sequence with N replaced by G, as qsidd reports it.

        This is what qsidd prints in its ``Base`` column, and what is written
        to ``input.analysed_sequence``.
        """
        return self.analysed.replace("N", "G")

    @property
    def sha256(self) -> str:
        """SHA-256 of :attr:`bases`, as 64 lowercase hex characters.

        Written to ``input.sequence_sha256``, for comparing which sequence two
        results came from.
        """
        return hashlib.sha256(self.bases.encode("ascii")).hexdigest()

    @property
    def n_positions(self) -> list[int]:
        """1-based positions of every N in the analysed sequence."""
        return [index + 1 for index, base in enumerate(self.analysed) if base == "N"]
