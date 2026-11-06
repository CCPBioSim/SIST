"""Validate structured SIST results.

A result is checked against the JSON schema, which covers its shape: types,
ranges, required keys and allowed names.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass
from importlib import resources
from typing import Any

from jsonschema import Draft202012Validator
from jsonschema.exceptions import ValidationError


@dataclass(frozen=True)
class ValidationIssue:
    """A single problem found in a SIST result.

    Attributes:
        path: Where the problem is, written as keys joined by ``/``, for example
            ``"profile/probability/melting"``. Empty for the top level.
        message: What is wrong.
    """

    path: str
    message: str

    def __str__(self) -> str:
        """Format the issue as ``path: message``."""
        return f"{self.path or '(top level)'}: {self.message}"


class ResultValidationError(ValueError):
    """Error for a SIST result that fails validation.

    Attributes:
        issues: Every problem found in the result.
    """

    def __init__(self, issues: Sequence[ValidationIssue]) -> None:
        """Build the error message from every problem found.

        Args:
            issues: Problems found in the result.
        """
        self.issues = list(issues)

        lines = "\n".join(f"  {issue}" for issue in self.issues)
        super().__init__(f"{len(self.issues)} problem(s) in SIST result:\n{lines}")


class ResultValidator:
    """Check SIST results against the result schema.

    Attributes:
        SCHEMA_FILE: File name of the schema inside ``sist/schemas``.
        schema: The result schema, loaded from the package.
    """

    SCHEMA_FILE = "sist-result-v1.schema.json"

    def __init__(self) -> None:
        """Load the schema and set up the JSON Schema validator."""
        self.schema = self.load_schema()

        self._validator = Draft202012Validator(self.schema)

    @classmethod
    def load_schema(cls) -> dict[str, Any]:
        """Load the result schema that ships with the package.

        Returns:
            The schema as a dictionary.
        """
        schema_text = (
            resources.files("sist")
            .joinpath("schemas", cls.SCHEMA_FILE)
            .read_text(encoding="utf-8")
        )
        schema: dict[str, Any] = json.loads(schema_text)

        return schema

    def validate(self, result: dict[str, Any]) -> list[ValidationIssue]:
        """Find every problem in a result.

        Args:
            result: A parsed SIST result.

        Returns:
            Every problem found. An empty list means the result is valid.
        """
        return [
            self._schema_issue(error) for error in self._validator.iter_errors(result)
        ]

    def check(self, result: dict[str, Any]) -> None:
        """Make sure a result is valid.

        Use this instead of :meth:`validate` when a bad result should stop the
        program.

        Args:
            result: A parsed SIST result.

        Raises:
            ResultValidationError: If the result has any problems.
        """
        issues = self.validate(result)

        if issues:
            raise ResultValidationError(issues)

    @staticmethod
    def _schema_issue(error: ValidationError) -> ValidationIssue:
        """Turn a schema error into an issue that names the real problem.

        Sections that can be null, such as ``profile``, are written in the schema
        as "one of these options". When such a section is wrong, the schema error
        only says that no option matched. Following the error down to the
        deepest point in the result nearly always reaches the actual bad value.

        Args:
            error: An error from the JSON Schema validator.

        Returns:
            The issue for the deepest underlying error.
        """
        while error.context:
            error = max(error.context, key=lambda inner: len(inner.absolute_path))

        path = "/".join(str(part) for part in error.absolute_path)

        return ValidationIssue(path, error.message)
