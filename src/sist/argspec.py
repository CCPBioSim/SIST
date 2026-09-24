"""Command-line argument specification for sist.

This module provides a declarative argument specification (`ARG_SPECS`)
used to build an `argparse.ArgumentParser` for the sist CLI's short and
long-form flags.

"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Any

ALGORITHM_CHOICES = ("M", "Z", "C", "A")


@dataclass(frozen=True)
class ArgSpec:
    """Argument specification used to build an argparse parser.

    Attributes:
        flags: Option strings, e.g. ``("-f", "--file")``.
        dest: Attribute name the parsed value is stored under.
        help: Help text shown in CLI usage.
        default: Default value if not provided via CLI.
        type: Python type for parsing, such as int, float, or str.
        action: Optional argparse action, such as "store_true".
        choices: Optional set of allowed values.
        required: Whether the argument must be provided.
    """

    dest: str
    help: str
    flags: tuple[str, ...]
    default: Any = None
    type: Any = None
    action: str | None = None
    choices: tuple[str, ...] | None = None
    required: bool = False


ARG_SPECS: tuple[ArgSpec, ...] = (
    ArgSpec(
        dest="file",
        flags=("-f", "--file"),
        type=str,
        required=True,
        help="Sequence file to analyse.",
    ),
    ArgSpec(
        dest="algorithm",
        flags=("-a", "--algorithm"),
        type=str,
        required=True,
        choices=ALGORITHM_CHOICES,
        help=(
            "Algorithm to run: M (melting), Z (Z-DNA), C (cruciform), or A "
            "(competition between all three). C and A require an Inverted "
            "Repeat Finder (IRF) executable on PATH."
        ),
    ),
    ArgSpec(
        dest="temperature",
        flags=("-T", "--temperature"),
        type=float,
        default=310.0,
        help="Temperature in Kelvin. Defaults to %(default)s.",
    ),
    ArgSpec(
        dest="superhelical_density",
        flags=("-s", "--superhelical-density"),
        type=float,
        default=0.06,
        help="Superhelical density. Defaults to %(default)s.",
    ),
    ArgSpec(
        dest="salt",
        flags=("-i", "--salt"),
        type=float,
        default=0.01,
        help="Ionic strength (salt concentration) in mol/L. Defaults to %(default)s.",
    ),
    ArgSpec(
        dest="energy_threshold",
        flags=("-th", "--energy-threshold"),
        type=float,
        default=12.0,
        help=(
            "Energy threshold for the transition search. Defaults to "
            "%(default)s; 10 is recommended for the competition algorithm."
        ),
    ),
    ArgSpec(
        dest="circular",
        flags=("-c", "--circular"),
        action="store_true",
        help="Treat the molecule as circular. Defaults to linear.",
    ),
    ArgSpec(
        dest="nearest_neighbor",
        flags=("-n", "--nearest-neighbor"),
        action="store_true",
        help=(
            "Use nearest-neighbor melting energetics instead of "
            "copolymeric. Ignored for the Z-DNA and cruciform algorithms."
        ),
    ),
    ArgSpec(
        dest="print_base_pair",
        flags=("-b", "--print-base-pair"),
        action="store_true",
        help="Print the base pair for each position.",
    ),
    ArgSpec(
        dest="print_parameters",
        flags=("-p", "--print-parameters"),
        action="store_true",
        help="Print algorithm parameters.",
    ),
    ArgSpec(
        dest="print_ensemble_average",
        flags=("-r", "--print-ensemble-average"),
        action="store_true",
        help="Print ensemble average results.",
    ),
    ArgSpec(
        dest="output_file",
        flags=("-o", "--output-file"),
        type=str,
        default=None,
        help="Write output to this file instead of printing to stdout.",
    ),
)


class SistArgumentParser:
    """Builds and validates the sist CLI's argument parser."""

    def __init__(self, arg_specs: tuple[ArgSpec, ...] | None = None) -> None:
        """Initialise the parser builder.

        Args:
            arg_specs: Optional override for argument specs. If omitted,
                uses `ARG_SPECS`.
        """
        self._arg_specs = arg_specs if arg_specs is not None else ARG_SPECS

    def build_parser(self) -> argparse.ArgumentParser:
        """Build an ArgumentParser from the argument specs.

        Returns:
            Configured argparse.ArgumentParser.
        """
        parser = argparse.ArgumentParser(
            prog="sist",
            description=(
                "SIST: stress-induced structural transition probabilities in "
                "superhelical DNA."
            ),
        )

        for spec in self._arg_specs:
            kwargs: dict[str, Any] = {"dest": spec.dest, "help": spec.help}

            if spec.action is not None:
                kwargs["action"] = spec.action
            else:
                kwargs["type"] = spec.type
                kwargs["default"] = spec.default

                if spec.choices is not None:
                    kwargs["choices"] = spec.choices

                if spec.required:
                    kwargs["required"] = True

            parser.add_argument(*spec.flags, **kwargs)

        return parser

    @staticmethod
    def validate(args: argparse.Namespace) -> None:
        """Validate parsed arguments against sensible runtime constraints.

        Args:
            args: Parsed CLI arguments.

        Raises:
            ValueError: If a parameter is invalid.
        """
        if args.temperature <= 0:
            raise ValueError(
                f"Invalid 'temperature': {args.temperature}. Temperature must "
                "be greater than 0 K."
            )

        if args.salt <= 0:
            raise ValueError(
                f"Invalid 'salt': {args.salt}. Salt must be greater than 0 mol/L."
            )

        if args.energy_threshold <= 0:
            raise ValueError(
                f"Invalid 'energy_threshold': {args.energy_threshold}. It must "
                "be greater than 0."
            )
