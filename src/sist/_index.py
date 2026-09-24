"""Shared helper for indexing with negative wraparound and no bounds errors.

Used by modules that need to look at neighbouring positions in an
alignment without special-casing the edges of the sequence.
"""

from __future__ import annotations

from collections.abc import Sequence


def index_or_none[T](sequence: Sequence[T], index: int) -> T | None:
    """Index a sequence, returning ``None`` instead of raising if out of range.

    Negative indices wrap from the end, as with normal Python indexing.

    Args:
        sequence: Sequence to index.
        index: Index to read.

    Returns:
        The element at ``index``, or ``None`` if out of range.
    """
    if -len(sequence) <= index < len(sequence):
        return sequence[index]

    return None
