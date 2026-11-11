"""Unit tests for sist.argspec."""

from __future__ import annotations

import argparse

import pytest

from sist.argspec import SistArgumentParser
from sist.transitions import ALGORITHMS


def _parse(*argv: str) -> argparse.Namespace:
    """Parse a command line with the sist parser.

    Args:
        *argv: The arguments, without the program name.

    Returns:
        The parsed arguments.
    """
    return SistArgumentParser().build_parser().parse_args(list(argv))


def test_build_parser_reads_flags_and_defaults() -> None:
    args = _parse("-f", "seq.fa", "-a", "A", "-c", "-th", "10")

    assert args.file == "seq.fa"
    assert args.algorithm == "A"
    assert args.circular is True
    assert args.energy_threshold == 10.0
    assert args.temperature == 310.0


def test_build_parser_offers_every_algorithm(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit):
        _parse("-f", "seq.fa", "-a", "X")

    error = capsys.readouterr().err

    assert all(flag in error for flag in ALGORITHMS)


@pytest.mark.parametrize(
    ("flag", "name"),
    [("-T", "temperature"), ("-i", "salt"), ("-th", "energy_threshold")],
)
def test_validate_rejects_values_that_are_not_positive(flag: str, name: str) -> None:
    args = _parse("-f", "seq.fa", "-a", "M", flag, "0")

    with pytest.raises(ValueError, match=name):
        SistArgumentParser.validate(args)
