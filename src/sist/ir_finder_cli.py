"""Deprecated console-script entry point for the old `IR_finder.pl` interface.

`IR_finder.pl` took positional arguments (``<temperature> <shape>
<sequence_file>``) rather than flags, and printed the resulting IR string to
stdout. This module reproduces that exact interface on top of
:class:`sist.ir_finder.IRFinder` for scripted workflows that still invoke it
by that name (see ``[project.scripts]`` in ``pyproject.toml``). It will be
removed in the next release; use the ``sist`` command instead.
"""

from __future__ import annotations

import logging
import sys

from sist.ir_finder import IRFinder

logger = logging.getLogger(__name__)

USAGE = "usage: IR_finder.pl <temperature> <shape> <sequence_file>"


def main(argv: list[str] | None = None) -> None:
    """Parse `IR_finder.pl`-style positional arguments and print the IR string.

    Args:
        argv: Argument list to parse, or None to use `sys.argv[1:]`.

    Raises:
        SystemExit: Exits with status code 1 if the arguments are missing,
            or on any unhandled exception.
    """
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    print(
        "WARNING: 'IR_finder.pl' is a deprecated alias for 'sist' and will "
        "be removed in the next release. Use 'sist' instead.",
        file=sys.stderr,
    )

    args = sys.argv[1:] if argv is None else argv

    if len(args) < 3:
        print(USAGE, file=sys.stderr)
        raise SystemExit(1)

    temperature, shape, sequence_file = args[0], args[1], args[2]

    try:
        finder = IRFinder(temperature=float(temperature), shape=shape)
        print(finder.compute_cruciform_energy_string(sequence_file))
    except Exception:
        logger.exception("Fatal error during IR_finder calculation")
        raise SystemExit(1) from None
