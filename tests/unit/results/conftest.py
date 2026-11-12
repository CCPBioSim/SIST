"""Fixtures for the sist.results tests."""

from __future__ import annotations

import pytest
from results_helpers import read_output

from sist.results.qsidd_output import QsiddOutput


@pytest.fixture
def melting_output() -> QsiddOutput:
    return read_output("melting.txt")


@pytest.fixture
def competition_output() -> QsiddOutput:
    return read_output("competition.txt")
