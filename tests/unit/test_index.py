"""Unit tests for sist._index."""

from __future__ import annotations

import pytest

from sist._index import index_or_none


@pytest.mark.parametrize(
    ("index", "expected"), [(0, "A"), (-1, "C"), (3, None), (-4, None)]
)
def test_index_or_none(index: int, expected: str | None) -> None:
    assert index_or_none("ABC", index) == expected
