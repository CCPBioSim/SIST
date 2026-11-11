"""Check a result against the schema and write it as JSON."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from sist.results.model import Result
from sist.schema import ResultValidator


class ResultWriter:
    """Check results against the schema before writing them.

    Attributes:
        validator: Validator every result must pass.
    """

    def __init__(self, validator: ResultValidator | None = None) -> None:
        """Set up the writer.

        Args:
            validator: Validator to use, or None for a new one.
        """
        self.validator = validator or ResultValidator()

    def check(self, result: Result) -> dict[str, Any]:
        """Turn a result into its schema shape and make sure it is valid.

        Args:
            result: The result to check.

        Returns:
            The result as a dictionary.

        Raises:
            ResultValidationError: If the result is not valid.
        """
        data = result.to_dict()
        self.validator.check(data)

        return data

    def write(self, result: Result, path: Path) -> None:
        """Check a result and write it as JSON.

        Args:
            result: The result to write.
            path: File to write it to.

        Raises:
            ResultValidationError: If the result is not valid.
        """
        path.write_text(self.to_json(self.check(result)), encoding="utf-8")

    @staticmethod
    def to_json(data: dict[str, Any]) -> str:
        """Write a result dictionary as JSON.

        Objects are indented, and arrays of plain values are written on one
        line.

        Args:
            data: The result as a dictionary.

        Returns:
            The JSON text, ending with a newline.

        Raises:
            ValueError: If the result holds NaN or Infinity.
        """
        text = json.dumps(data, indent=2, allow_nan=False, ensure_ascii=False)
        value = r'(?:-?[0-9][0-9.eE+-]*|null|true|false|"[^"\n]*")'
        flat_array = re.compile(rf"\[\n(\s+{value}(?:,\n\s+{value})*)\n\s*\]")

        def join(match: re.Match[str]) -> str:
            items = [item.strip() for item in match.group(1).split(",\n")]
            return "[" + ", ".join(items) + "]"

        return flat_array.sub(join, text) + "\n"
