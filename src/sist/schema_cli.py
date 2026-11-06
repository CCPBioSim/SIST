"""Command-line entry point for validating SIST result files.

Run it on one or more files to check each one and list any problems::

    sist-validate result.json [result.json ...]

It is installed as the ``sist-validate`` console script (see
``[project.scripts]`` in ``pyproject.toml``).
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Any

from sist.schema import ResultValidator

logger = logging.getLogger(__name__)

USAGE = "usage: sist-validate RESULT.json [RESULT.json ...]"


def main(argv: list[str] | None = None) -> None:
    """Check each SIST result file and print what was found.

    Every file is checked before exiting, so one bad file does not hide
    problems in the others.

    Args:
        argv: Result file paths, or None to use ``sys.argv[1:]``.

    Raises:
        SystemExit: With status code 1 if any file cannot be read or is not
            valid, or on any unhandled exception, and 2 if no files are given.
    """
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

    paths = sys.argv[1:] if argv is None else argv

    if not paths:
        print(USAGE, file=sys.stderr)
        raise SystemExit(2)

    all_valid = True

    try:
        validator = ResultValidator()

        for path in paths:
            if not _report(validator, Path(path)):
                all_valid = False
    except Exception:
        logger.exception("Fatal error while validating SIST results")
        raise SystemExit(1) from None

    if not all_valid:
        raise SystemExit(1)


def _report(validator: ResultValidator, path: Path) -> bool:
    """Check one result file and print the outcome.

    Args:
        validator: Validator to check the result with.
        path: Path to the result file.

    Returns:
        True if the file could be read and is valid, otherwise False.
    """
    try:
        result = _read_result(path)
    except (OSError, ValueError) as error:
        print(f"{path}: cannot read: {error}")
        return False

    issues = validator.validate(result)

    if not issues:
        print(f"{path}: valid")
        return True

    print(f"{path}: {len(issues)} problem(s)")

    for issue in issues:
        print(f"  {issue}")

    return False


def _read_result(path: Path) -> dict[str, Any]:
    """Read a result file, refusing anything that is not strict JSON.

    Python's json module accepts ``NaN`` and ``Infinity`` by default, even
    though they are not valid JSON and other tools cannot read them. Refusing
    them means a file that passes here can be read anywhere.

    Args:
        path: Path to the result file.

    Returns:
        The parsed result.

    Raises:
        ValueError: If the file is not valid JSON.
    """
    result: dict[str, Any] = json.loads(
        path.read_text(encoding="utf-8"), parse_constant=_reject_constant
    )

    return result


def _reject_constant(constant: str) -> None:
    """Refuse a value that is not valid JSON, such as ``NaN``.

    Args:
        constant: The value found.

    Raises:
        ValueError: Always.
    """
    raise ValueError(f"{constant} is not valid JSON")
