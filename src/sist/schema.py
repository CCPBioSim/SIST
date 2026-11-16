"""Validate structured SIST results.

A result is checked in two passes:

1. The JSON schema, ``sist-result-v1.schema.json``, checks its shape: types,
   ranges, required keys and allowed names.
2. The cross-checks compare one part of the result with another: the
   transitions against the algorithm, each per-position array against the
   sequence length, the state counts against their total, and the partition
   function against its log.

Each problem is reported as a :class:`ValidationIssue`, which prints as its
path and message::

    for issue in ResultValidator().validate(result):
        print(issue)

For a melting profile one entry too long, this prints a line starting
``transitions/melting/profile/probability:``.
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
        path: Where the problem is, written as keys and list indexes joined by
            ``/``, for example ``"transitions/melting/profile/probability/2"``.
            Empty for the top level.
        message: What is wrong, such as ``"1.5 is greater than the maximum of
            1"``.
    """

    path: str
    message: str

    def __str__(self) -> str:
        """Format the issue as ``path: message``, or ``(top level): message``."""
        return f"{self.path or '(top level)'}: {self.message}"


class ResultValidationError(ValueError):
    """Error for a SIST result that fails validation.

    The message lists every problem, one per line::

        2 problem(s) in SIST result:
          schema_version: '2.0.0' does not match ...
          input/analysed_sequence: has 6 letters, expected 5

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

    :meth:`validate` returns every problem found, and :meth:`check` raises
    :class:`ResultValidationError` if there are any. One validator can check
    any number of results.

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

        The schema is read from ``sist/schemas/sist-result-v1.schema.json``
        inside the installed package.

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

        The result is checked against the schema first. If the schema passes,
        the cross-checks run as well. A result that fails the schema reports
        only the schema problems.

        Args:
            result: A parsed SIST result, such as the output of ``json.load``.

        Returns:
            Every problem found, schema problems first. An empty list means the
            result is valid.
        """
        issues = [
            self._schema_issue(error) for error in self._validator.iter_errors(result)
        ]

        if issues:
            return issues

        return self._cross_check(result)

    def check(self, result: dict[str, Any]) -> None:
        """Raise an error listing every problem in a result.

        Runs :meth:`validate`, and does nothing if it finds no problems.

        Args:
            result: A parsed SIST result.

        Raises:
            ResultValidationError: If the result has any problems. Its
                ``issues`` attribute holds them all.
        """
        issues = self.validate(result)

        if issues:
            raise ResultValidationError(issues)

    @staticmethod
    def _schema_issue(error: ValidationError) -> ValidationIssue:
        """Turn a schema error into an issue for its deepest underlying error.

        For a section with more than one allowed shape, such as ``transitions``
        (results or null), the schema error holds one inner error per shape.
        This follows the inner errors down to the one with the longest path in
        the result, so a bad probability inside ``transitions`` is reported at
        ``transitions/melting/profile/probability/2``.

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

        Every result gets the timestamp and transition checks. A successful
        result also gets the pairwise, length, state count and partition
        function checks, which read its ``summary`` and ``transitions``.

        Args:
            result: A result that has passed the schema.

        Returns:
            Every problem found, in the order the checks run.
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

        Reads ``provenance.started_at`` and ``provenance.finished_at``, and
        parses each with :meth:`datetime.fromisoformat`. A null time is
        accepted.

        Args:
            result: A result that has passed the schema.

        Returns:
            An issue at ``provenance/started_at`` or ``provenance/finished_at``
            for each time that does not parse.
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

        Two comparisons are made:

        * ``calculation.transitions`` against the transitions the algorithm in
          ``calculation.algorithm`` calculates, from :data:`ALGORITHMS`.
        * The keys of ``transitions`` against ``calculation.transitions``. This
          is skipped for a failed run, where ``transitions`` is null.

        Args:
            result: A result that has passed the schema.

        Returns:
            An issue at ``calculation/transitions`` if they do not match the
            algorithm, and issues at ``transitions`` naming any extra or
            missing transitions.
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

        Reads ``summary.pairwise``, if present. Each of its three matrices is
        read as ``matrix[first][second]``, and every pair is checked against
        ``calculation.transitions``.

        Args:
            result: A successful result that has passed the schema.

        Returns:
            One issue at ``summary/pairwise`` if the calculation has only one
            transition. Otherwise an issue at
            ``summary/pairwise/<matrix>/<first>/<second>`` for each pair that
            repeats a transition or names one not in the calculation.
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

        Compares ``input.sequence_length_bp`` with the length of
        ``input.analysed_sequence`` and of every array in each transition's
        ``profile``.

        Args:
            result: A successful result that has passed the schema.

        Returns:
            One issue at ``input/sequence_length_bp`` if the length is null.
            Otherwise an issue for each value of the wrong length, giving its
            length and the expected one.
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

        Sums ``summary.state_counts.by_number_of_runs`` and compares the sum
        with ``summary.state_counts.total``.

        Args:
            result: A successful result that has passed the schema.

        Returns:
            An issue at ``summary/state_counts`` giving both numbers if they
            differ, otherwise an empty list.
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

        A null ``summary.partition_function`` is accepted. Otherwise its log is
        compared with ``summary.ln_partition_function``, to within
        :attr:`LN_PARTITION_FUNCTION_RTOL`.

        Args:
            result: A successful result that has passed the schema.

        Returns:
            An issue at ``summary/partition_function`` if it is 0 or below, or
            at ``summary/ln_partition_function`` if the two do not match.
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
