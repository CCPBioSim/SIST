"""Validate structured SIST results.

A result is checked in two passes. The JSON schema checks its shape: types,
ranges, required keys and allowed names. The cross-checks then compare one part
of the result with another, such as each per-position array against the
sequence length.
"""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from importlib import resources
from typing import Any

from jsonschema import Draft202012Validator
from jsonschema.exceptions import ValidationError

from sist.transitions import ALGORITHMS


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
    """Check SIST results against the result schema and the cross-checks.

    Attributes:
        SCHEMA_FILE: File name of the schema inside ``sist/schemas``.
        LN_PARTITION_FUNCTION_RTOL: Relative tolerance when comparing
            ``ln_partition_function`` with the log of ``partition_function``.
        schema: The result schema, loaded from the package.
    """

    SCHEMA_FILE = "sist-result-v1.schema.json"

    LN_PARTITION_FUNCTION_RTOL = 1e-9

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

        The cross-checks run only when the schema passes. A result that fails
        the schema reports only the schema problems.

        Args:
            result: A parsed SIST result.

        Returns:
            Every problem found. An empty list means the result is valid.
        """
        issues = [
            self._schema_issue(error) for error in self._validator.iter_errors(result)
        ]

        if issues:
            return issues

        return self._cross_check(result)

    def check(self, result: dict[str, Any]) -> None:
        """Raise an error listing every problem in a result.

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
        """Turn a schema error into an issue for its deepest underlying error.

        For a section with more than one allowed shape, such as ``transitions``
        (results or null), this follows the error down to the deepest point in
        the result.

        Args:
            error: An error from the JSON Schema validator.

        Returns:
            The issue for the deepest underlying error.
        """
        while error.context:
            error = max(error.context, key=lambda inner: len(inner.absolute_path))

        path = "/".join(str(part) for part in error.absolute_path)

        return ValidationIssue(path, error.message)

    def _cross_check(self, result: dict[str, Any]) -> list[ValidationIssue]:
        """Run the checks that compare one part of a result with another.

        A failed run gets the timestamp and transition checks. A successful run
        gets every check.

        Args:
            result: A result that has passed the schema.

        Returns:
            Every problem found.
        """
        issues = []

        issues.extend(self._check_timestamps(result))
        issues.extend(self._check_transitions(result))

        if result["status"] == "success":
            issues.extend(self._check_pairwise(result))
            issues.extend(self._check_lengths(result))
            issues.extend(self._check_state_counts(result))
            issues.extend(self._check_partition_function(result))

        return issues

    @staticmethod
    def _check_timestamps(result: dict[str, Any]) -> list[ValidationIssue]:
        """Check that the start and finish times are RFC 3339 date-times.

        Args:
            result: A result that has passed the schema.

        Returns:
            A problem for each time that is not an RFC 3339 date-time.
        """
        issues = []
        provenance = result["provenance"]

        for key in ("started_at", "finished_at"):
            timestamp = provenance[key]

            if timestamp is None:
                continue

            try:
                datetime.fromisoformat(timestamp)
            except ValueError:
                issues.append(
                    ValidationIssue(
                        f"provenance/{key}",
                        f"{timestamp!r} is not an RFC 3339 date-time",
                    )
                )

        return issues

    @staticmethod
    def _check_transitions(result: dict[str, Any]) -> list[ValidationIssue]:
        """Check the transitions against the algorithm and the results.

        ``calculation.transitions`` is compared with what the algorithm
        calculates, and ``transitions`` with ``calculation.transitions``.

        Args:
            result: A result that has passed the schema.

        Returns:
            A problem if the transitions do not match the algorithm, and one if
            the results have extra or missing transitions.
        """
        issues = []
        calculation = result["calculation"]
        transitions = calculation["transitions"]
        algorithm = ALGORITHMS.get(calculation["algorithm"])

        if algorithm is not None and tuple(transitions) != algorithm.transition_names:
            issues.append(
                ValidationIssue(
                    "calculation/transitions",
                    f"algorithm {algorithm.flag} calculates "
                    f"{list(algorithm.transition_names)}, not {transitions}",
                )
            )

        results = result["transitions"]

        if results is None:
            return issues

        extra = sorted(set(results) - set(transitions))
        missing = sorted(set(transitions) - set(results))

        if extra:
            issues.append(
                ValidationIssue(
                    "transitions",
                    f"has transitions not in calculation.transitions: {extra}",
                )
            )

        if missing:
            issues.append(
                ValidationIssue("transitions", f"is missing transitions: {missing}")
            )

        return issues

    @staticmethod
    def _check_pairwise(result: dict[str, Any]) -> list[ValidationIssue]:
        """Check that pairwise probabilities only pair different, known transitions.

        Args:
            result: A successful result that has passed the schema.

        Returns:
            A problem if the calculation has only one transition, and one for
            each pair that repeats a transition or uses one not in the
            calculation.
        """
        pairwise = result["summary"].get("pairwise")

        if pairwise is None:
            return []

        transitions = set(result["calculation"]["transitions"])

        if len(transitions) < 2:
            return [
                ValidationIssue("summary/pairwise", "needs more than one transition")
            ]

        issues = []

        for name, matrix in pairwise.items():
            for first, row in matrix.items():
                for second in row:
                    path = f"summary/pairwise/{name}/{first}/{second}"

                    if first == second:
                        issues.append(
                            ValidationIssue(path, "pairs a transition with itself")
                        )
                    elif {first, second} - transitions:
                        issues.append(
                            ValidationIssue(
                                path,
                                "uses a transition not in calculation.transitions",
                            )
                        )

        return issues

    @staticmethod
    def _check_lengths(result: dict[str, Any]) -> list[ValidationIssue]:
        """Check that every per-position value has one entry per base pair.

        Args:
            result: A successful result that has passed the schema.

        Returns:
            A problem for each array, or the analysed sequence, of the wrong
            length, or a single problem if the sequence length is missing.
        """
        length = result["input"]["sequence_length_bp"]

        if length is None:
            return [
                ValidationIssue(
                    "input/sequence_length_bp", "must be set for a successful run"
                )
            ]

        issues = []
        analysed_sequence = result["input"].get("analysed_sequence")

        if analysed_sequence is not None and len(analysed_sequence) != length:
            issues.append(
                ValidationIssue(
                    "input/analysed_sequence",
                    f"has {len(analysed_sequence)} letters, expected {length}",
                )
            )

        for transition, transition_result in result["transitions"].items():
            for name, values in transition_result["profile"].items():
                if len(values) != length:
                    issues.append(
                        ValidationIssue(
                            f"transitions/{transition}/profile/{name}",
                            f"has {len(values)} entries, expected {length}",
                        )
                    )

        return issues

    @staticmethod
    def _check_state_counts(result: dict[str, Any]) -> list[ValidationIssue]:
        """Check that the state counts add up to the total.

        Args:
            result: A successful result that has passed the schema.

        Returns:
            A problem if the counts do not add up to the total.
        """
        state_counts = result["summary"]["state_counts"]
        counted = sum(state_counts["by_number_of_runs"])

        if counted == state_counts["total"]:
            return []

        return [
            ValidationIssue(
                "summary/state_counts",
                f"by_number_of_runs adds up to {counted}, "
                f"but total is {state_counts['total']}",
            )
        ]

    @classmethod
    def _check_partition_function(cls, result: dict[str, Any]) -> list[ValidationIssue]:
        """Check that the partition function is above 0 and matches its natural log.

        Args:
            result: A successful result that has passed the schema.

        Returns:
            A problem if the partition function is not above 0, or does not
            match ``ln_partition_function``.
        """
        summary = result["summary"]
        partition_function = summary["partition_function"]

        if partition_function is None:
            return []

        if partition_function <= 0:
            return [
                ValidationIssue(
                    "summary/partition_function",
                    "must be above 0; use null if it is too small to store",
                )
            ]

        expected = math.log(partition_function)

        if math.isclose(
            expected,
            summary["ln_partition_function"],
            rel_tol=cls.LN_PARTITION_FUNCTION_RTOL,
        ):
            return []

        return [
            ValidationIssue(
                "summary/ln_partition_function",
                f"does not match log(partition_function) = {expected}",
            )
        ]
