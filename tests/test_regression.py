from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest

from conftest import CompetitionRun


pytestmark = pytest.mark.regression


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

REFERENCE_OUTPUT = (
    REPOSITORY_ROOT
    / "tests"
    / "reference"
    / "v1.0.0"
    / "competition.rebuilt.txt"
)

# The baseline reproduced all shared printed values exactly.
RELATIVE_TOLERANCE = 0.0
ABSOLUTE_TOLERANCE = 0.0


@dataclass(frozen=True)
class ProfileRow:
    """One position from the SIST probability profile."""

    base: str
    p_melt: float
    p_z: float
    p_cruciform: float


@dataclass(frozen=True)
class ParsedOutput:
    """Parsed deterministic content from a SIST result."""

    metadata: tuple[str, ...]
    metrics: dict[str, float]
    profile: dict[int, ProfileRow]


def parse_numeric_value(value: str) -> float | None:
    """
    Parse the first token as a number.

    This handles values followed by units, such as:
    '-0.06'
    '310 K'
    '1.4748e-54'
    """

    first_token = value.strip().split()[0]

    try:
        return float(first_token)
    except ValueError:
        return None


def parse_output(path: Path) -> ParsedOutput:
    """Parse deterministic metrics and profile rows from SIST output."""

    metadata: list[str] = []
    metrics: dict[str, float] = {}
    profile: dict[int, ProfileRow] = {}

    reading_profile = False

    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        line = raw_line.strip()

        if not line:
            continue

        # Runtime varies between executions and is not scientific output.
        if line.startswith("Run time ="):
            continue

        if line.startswith("Position"):
            reading_profile = True
            continue

        if reading_profile:
            columns = line.split()

            if len(columns) != 5:
                raise AssertionError(
                    f"{path}:{line_number}: expected five profile columns, "
                    f"found {len(columns)}: {line!r}"
                )

            position_text, base, p_melt, p_z, p_cruciform = columns

            position = int(position_text)

            profile[position] = ProfileRow(
                base=base,
                p_melt=float(p_melt),
                p_z=float(p_z),
                p_cruciform=float(p_cruciform),
            )

            continue

        if line.startswith("Scaling factor "):
            value = parse_numeric_value(
                line.removeprefix("Scaling factor ")
            )

            if value is None:
                raise AssertionError(
                    f"{path}:{line_number}: could not parse scaling factor"
                )

            metrics["Scaling factor"] = value
            continue

        # Some lines contain more than one key/value pair separated by ';'.
        parsed_segment = False

        for segment in line.split(";"):
            if " = " not in segment:
                continue

            key, raw_value = segment.split(" = ", maxsplit=1)
            numeric_value = parse_numeric_value(raw_value)

            if numeric_value is None:
                continue

            metrics[key.strip()] = numeric_value
            parsed_segment = True

        if not parsed_segment:
            metadata.append(line)

    if not profile:
        raise AssertionError(f"No position profile was found in {path}")

    return ParsedOutput(
        metadata=tuple(metadata),
        metrics=metrics,
        profile=profile,
    )


def assert_number_matches(
    *,
    name: str,
    expected: float,
    actual: float,
) -> None:
    """Compare a numerical value with the documented baseline tolerance."""

    assert actual == pytest.approx(
        expected,
        rel=RELATIVE_TOLERANCE,
        abs=ABSOLUTE_TOLERANCE,
    ), (
        f"{name} changed: expected {expected!r}, "
        f"actual {actual!r}"
    )


def test_competition_scientific_results_match_baseline(
    competition_run: CompetitionRun,
) -> None:
    """The rebuilt calculation should reproduce the baseline results."""

    expected = parse_output(REFERENCE_OUTPUT)
    actual = parse_output(competition_run.output_path)

    assert actual.metadata == expected.metadata, (
        "Non-numerical output metadata changed:\n"
        f"expected: {expected.metadata!r}\n"
        f"actual:   {actual.metadata!r}"
    )

    assert actual.metrics.keys() == expected.metrics.keys(), (
        "The set of reported scientific metrics changed:\n"
        f"missing: {expected.metrics.keys() - actual.metrics.keys()}\n"
        f"extra:   {actual.metrics.keys() - expected.metrics.keys()}"
    )

    for metric_name, expected_value in expected.metrics.items():
        assert_number_matches(
            name=metric_name,
            expected=expected_value,
            actual=actual.metrics[metric_name],
        )

    assert actual.profile.keys() == expected.profile.keys(), (
        "The set of reported sequence positions changed"
    )

    for position, expected_row in expected.profile.items():
        actual_row = actual.profile[position]

        assert actual_row.base == expected_row.base, (
            f"Base changed at position {position}: "
            f"expected {expected_row.base!r}, "
            f"actual {actual_row.base!r}"
        )

        assert_number_matches(
            name=f"position {position} P_melt",
            expected=expected_row.p_melt,
            actual=actual_row.p_melt,
        )

        assert_number_matches(
            name=f"position {position} P_Z",
            expected=expected_row.p_z,
            actual=actual_row.p_z,
        )

        assert_number_matches(
            name=f"position {position} P_cruciform",
            expected=expected_row.p_cruciform,
            actual=actual_row.p_cruciform,
        )
