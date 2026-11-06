"""Validate structured SIST results.

A result is checked in two passes. The JSON schema checks its shape: types,
ranges, required keys and allowed names. Then a set of cross-checks compare one
part of the result with another, such as making sure every per-position array is
as long as the sequence. JSON Schema cannot express those comparisons, so they
are written in Python here.
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
        ALGORITHM_TRANSITIONS: Transitions calculated by each algorithm, in the
            order SIST reports them.
        COMPLETE_TRANSITION_SECTIONS: Sections that need an entry for every
            transition in the calculation, given as the keys to reach them.
        PARTIAL_TRANSITION_SECTIONS: Sections that may cover only some of the
            transitions. G(x) is one, because qsidd only reports it for melting.
        LN_PARTITION_FUNCTION_RTOL: Relative tolerance when comparing
            ``ln_partition_function`` with the log of ``partition_function``.
        schema: The result schema, loaded from the package.
    """

    SCHEMA_FILE = "sist-result-v1.schema.json"

    ALGORITHM_TRANSITIONS: dict[str, tuple[str, ...]] = {
        "M": ("melting",),
        "Z": ("z_dna",),
        "C": ("cruciform",),
        "A": ("melting", "z_dna", "cruciform"),
    }

    COMPLETE_TRANSITION_SECTIONS: tuple[tuple[str, ...], ...] = (
        ("summary", "transitions"),
        ("profile", "probability"),
        ("derived_parameters", "nucleation_energy_kcal_per_mol"),
    )

    PARTIAL_TRANSITION_SECTIONS: tuple[tuple[str, ...], ...] = (
        ("profile", "free_energy_kcal_per_mol"),
    )

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

        The cross-checks rely on the shape the schema guarantees, so they only
        run once the schema passes. A result that fails the schema reports
        just the schema problems.

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

    def _cross_check(self, result: dict[str, Any]) -> list[ValidationIssue]:
        """Compare parts of a result that the schema cannot check against each other.

        A failed run has no summary or profile, so only the checks that apply to
        every run are made for it.

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
        """Make sure the start and finish times are real dates.

        The schema marks them as dates, but the jsonschema package skips that
        check unless an optional extra is installed. Checking them here means
        they are always checked.

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

    def _check_transitions(self, result: dict[str, Any]) -> list[ValidationIssue]:
        """Make sure every section covers the transitions the calculation names.

        ``calculation.transitions`` is the single source of truth. It has to
        match the algorithm, and every section keyed by transition has to match
        it.

        Args:
            result: A result that has passed the schema.

        Returns:
            A problem if the transitions do not match the algorithm, and one for
            each section with extra or missing transitions.
        """
        issues = []
        calculation = result["calculation"]
        transitions = calculation["transitions"]
        expected = self.ALGORITHM_TRANSITIONS.get(calculation["algorithm"])

        if expected is not None and tuple(transitions) != expected:
            issues.append(
                ValidationIssue(
                    "calculation/transitions",
                    f"algorithm {calculation['algorithm']} calculates "
                    f"{list(expected)}, not {transitions}",
                )
            )

        sections = [(keys, True) for keys in self.COMPLETE_TRANSITION_SECTIONS]
        sections += [(keys, False) for keys in self.PARTIAL_TRANSITION_SECTIONS]
        wanted = set(transitions)

        for keys, must_cover_all in sections:
            section = self._get(result, *keys)

            if section is None:
                continue

            path = "/".join(keys)
            extra = sorted(set(section) - wanted)
            missing = sorted(wanted - set(section))

            if extra:
                issues.append(
                    ValidationIssue(
                        path, f"has transitions not in calculation.transitions: {extra}"
                    )
                )

            if must_cover_all and missing:
                issues.append(
                    ValidationIssue(path, f"is missing transitions: {missing}")
                )

        return issues

    @staticmethod
    def _check_pairwise(result: dict[str, Any]) -> list[ValidationIssue]:
        """Make sure pairwise probabilities only pair different, known transitions.

        The schema allows any transition name in either position, because it
        cannot compare a key with the key above it. So it cannot stop a
        transition being paired with itself, which would not be a valid
        probability.

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
        """Make sure every per-position value has one entry per base pair.

        Args:
            result: A successful result that has passed the schema.

        Returns:
            A problem for each array, or the base string, of the wrong length.
        """
        issues = []
        length = result["input"]["sequence_length_bp"]
        profile = result["profile"]

        if "base" in profile and len(profile["base"]) != length:
            issues.append(
                ValidationIssue(
                    "profile/base",
                    f"has {len(profile['base'])} letters, expected {length}",
                )
            )

        for group in ("probability", "free_energy_kcal_per_mol"):
            for transition, values in profile.get(group, {}).items():
                if len(values) != length:
                    issues.append(
                        ValidationIssue(
                            f"profile/{group}/{transition}",
                            f"has {len(values)} entries, expected {length}",
                        )
                    )

        return issues

    @staticmethod
    def _check_state_counts(result: dict[str, Any]) -> list[ValidationIssue]:
        """Make sure the state counts add up to the total.

        This catches counts that have lost digits. qsidd prints large counts in
        scientific notation, so 4908539 comes out as ``4.90854e+06``. That reads
        back as 4908540, a whole number, so it passes the schema's integer check.
        The total is always printed in full, so comparing against it shows the
        difference.

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
        """Make sure the partition function matches its natural log.

        A partition function too small to store should be written as null, not
        0, so a value of 0 or below is a problem.

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

    @staticmethod
    def _get(result: dict[str, Any], *keys: str) -> Any:
        """Look up a value nested inside a result.

        Args:
            result: A parsed SIST result.
            *keys: Keys to follow, outermost first.

        Returns:
            The value, or None if any key along the way is missing or null.
        """
        value: Any = result

        for key in keys:
            if not isinstance(value, dict):
                return None

            value = value.get(key)

        return value
