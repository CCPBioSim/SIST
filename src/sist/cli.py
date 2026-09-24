"""Command-line entry point for sist.

The entry point is intentionally small and only responsible for:
  1) Building the parser and parsing arguments.
  2) Validating those arguments.
  3) Constructing a SistRunner and running it.
  4) Handling fatal errors with a non-zero exit code.

It is also installed under the deprecated console-script alias ``master.pl``
(see ``[project.scripts]`` in ``pyproject.toml``) for scripted workflows that
still invoke it by that name. That alias will be removed in the next release.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from sist.argspec import SistArgumentParser
from sist.runner import SistRunner

logger = logging.getLogger(__name__)

DEPRECATED_ALIAS = "master.pl"


def main(argv: list[str] | None = None) -> None:
    """Parse CLI arguments and run a SIST calculation.

    Args:
        argv: Argument list to parse, or None to use `sys.argv[1:]`.

    Raises:
        SystemExit: Exits with status code 1 on any unhandled exception.
    """
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    if Path(sys.argv[0]).name == DEPRECATED_ALIAS:
        print(
            f"WARNING: '{DEPRECATED_ALIAS}' is a deprecated alias for 'sist' "
            "and will be removed in the next release. Use 'sist' instead.",
            file=sys.stderr,
        )

    argument_parser = SistArgumentParser()
    parser = argument_parser.build_parser()
    args = parser.parse_args(argv)

    try:
        argument_parser.validate(args)
        runner = SistRunner(
            file=args.file,
            algorithm=args.algorithm,
            temperature=args.temperature,
            superhelical_density=args.superhelical_density,
            salt=args.salt,
            energy_threshold=args.energy_threshold,
            circular=args.circular,
            nearest_neighbor=args.nearest_neighbor,
            print_base_pair=args.print_base_pair,
            print_parameters=args.print_parameters,
            print_ensemble_average=args.print_ensemble_average,
            output_file=args.output_file,
        )
        runner.run()
    except Exception:
        logger.exception("Fatal error during SIST calculation")
        raise SystemExit(1) from None
