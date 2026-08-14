from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest

from conftest import SistRun


pytestmark = pytest.mark.regression


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_VERSION = "v1.0.0"

REFERENCE_DIRECTORY = (
    REPOSITORY_ROOT
    / "tests"
    / "reference"
    / REFERENCE_VERSION
)

COMPETITION_REFERENCE = (
    REFERENCE_DIRECTORY
    / "competition.rebuilt.txt"
)

# The baselines reproduced all shared printed values exactly.
RELATIVE_TOLERANCE = 0.0
ABSOLUTE_TOLERANCE = 0.0


@dataclass(frozen=True)
class CompetitionProfileRow:
    """One position from the SIST competition probability profile."""

    base: str
    p_melt: float
    p_z: float
    p_cruciform: float


@dataclass(frozen=True)
class ParsedCompetitionOutput:
    """Parsed deterministic content from a SIST competition result."""

    metadata: tuple[str, ...]
    metrics: dict[str, float]
    profile: dict[int, CompetitionProfileRow]


@dataclass(frozen=True)
class TransitionProfileRow:
    """One position from an individual SIST transition profile."""

    base: str
    probability: float
    energy: float | None


@dataclass(frozen=True)
class ParsedTransitionOutput:
    """Parsed deterministic content from an individual SIST transition."""

    metadata: tuple[str, ...]
    metrics: dict[str, float]
    profile: dict[int, TransitionProfileRow]


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


def parse_competition_output(path: Path) -> ParsedCompetitionOutput:
    """Parse deterministic metrics and profile rows from competition output."""

    metadata: list[str] = []
    metrics: dict[str, float] = {}
    profile: dict[int, CompetitionProfileRow] = {}

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
            columns = line.split()

            expected_columns = (
                "Position",
                "Base",
                "P_melt",
                "P_Z",
                "P_cruciform",
            )

            if tuple(columns) != expected_columns:
                raise AssertionError(
                    f"{path}:{line_number}: unexpected competition "
                    f"profile header: {line!r}"
                )

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

            profile[position] = CompetitionProfileRow(
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

    return ParsedCompetitionOutput(
        metadata=tuple(metadata),
        metrics=metrics,
        profile=profile,
    )


def parse_transition_output(path: Path) -> ParsedTransitionOutput:
    """Parse deterministic metrics and profile rows from transition output."""

    metadata: list[str] = []
    metrics: dict[str, float] = {}
    profile: dict[int, TransitionProfileRow] = {}

    reading_profile = False
    profile_has_energy: bool | None = None

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
            columns = line.split()

            if columns == ["Position", "Base", "P(x)"]:
                profile_has_energy = False
            elif columns == ["Position", "Base", "P(x)", "G(x)"]:
                profile_has_energy = True
            else:
                raise AssertionError(
                    f"{path}:{line_number}: unexpected transition "
                    f"profile header: {line!r}"
                )

            reading_profile = True
            continue

        if reading_profile:
            if profile_has_energy is None:
                raise AssertionError(
                    f"{path}:{line_number}: profile format was not determined"
                )

            columns = line.split()
            expected_column_count = 4 if profile_has_energy else 3

            if len(columns) != expected_column_count:
                raise AssertionError(
                    f"{path}:{line_number}: expected "
                    f"{expected_column_count} profile columns, "
                    f"found {len(columns)}: {line!r}"
                )

            position = int(columns[0])
            base = columns[1]
            probability = float(columns[2])

            energy = (
                float(columns[3])
                if profile_has_energy
                else None
            )

            profile[position] = TransitionProfileRow(
                base=base,
                probability=probability,
                energy=energy,
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

    return ParsedTransitionOutput(
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


def assert_metadata_matches(
    *,
    name: str,
    expected: tuple[str, ...],
    actual: tuple[str, ...],
) -> None:
    """Compare deterministic non-numerical output metadata."""

    assert actual == expected, (
        f"{name} non-numerical metadata changed:\n"
        f"expected: {expected!r}\n"
        f"actual:   {actual!r}"
    )


def assert_metrics_match(
    *,
    name: str,
    expected: dict[str, float],
    actual: dict[str, float],
) -> None:
    """Compare all deterministic scientific metrics."""

    assert actual.keys() == expected.keys(), (
        f"{name} reported scientific metrics changed:\n"
        f"missing: {expected.keys() - actual.keys()}\n"
        f"extra:   {actual.keys() - expected.keys()}"
    )

    for metric_name, expected_value in expected.items():
        assert_number_matches(
            name=f"{name} {metric_name}",
            expected=expected_value,
            actual=actual[metric_name],
        )


def test_competition_scientific_results_match_baseline(
    competition_run: SistRun,
) -> None:
    """The competition calculation should reproduce its 1.0.0 baseline."""

    process = competition_run.process

    assert process.returncode == 0, (
        f"Competition calculation failed with exit code "
        f"{process.returncode}\n\n"
        f"stdout:\n{process.stdout}\n\n"
        f"stderr:\n{process.stderr}"
    )

    expected = parse_competition_output(COMPETITION_REFERENCE)
    actual = parse_competition_output(competition_run.output_path)

    assert_metadata_matches(
        name="competition",
        expected=expected.metadata,
        actual=actual.metadata,
    )

    assert_metrics_match(
        name="competition",
        expected=expected.metrics,
        actual=actual.metrics,
    )

    assert actual.profile.keys() == expected.profile.keys(), (
        "Competition reported sequence positions changed"
    )

    for position, expected_row in expected.profile.items():
        actual_row = actual.profile[position]

        assert actual_row.base == expected_row.base, (
            f"Competition base changed at position {position}: "
            f"expected {expected_row.base!r}, "
            f"actual {actual_row.base!r}"
        )

        assert_number_matches(
            name=f"competition position {position} P_melt",
            expected=expected_row.p_melt,
            actual=actual_row.p_melt,
        )

        assert_number_matches(
            name=f"competition position {position} P_Z",
            expected=expected_row.p_z,
            actual=actual_row.p_z,
        )

        assert_number_matches(
            name=f"competition position {position} P_cruciform",
            expected=expected_row.p_cruciform,
            actual=actual_row.p_cruciform,
        )


def test_transition_scientific_results_match_baseline(
    transition_run: SistRun,
) -> None:
    """Each individual transition should reproduce its 1.0.0 baseline."""

    process = transition_run.process

    assert process.returncode == 0, (
        f"{transition_run.name} calculation failed with exit code "
        f"{process.returncode}\n\n"
        f"stdout:\n{process.stdout}\n\n"
        f"stderr:\n{process.stderr}"
    )

    reference_output = (
        REFERENCE_DIRECTORY
        / f"{transition_run.name}.txt"
    )

    assert reference_output.is_file(), (
        f"Reference output does not exist: {reference_output}"
    )

    expected = parse_transition_output(reference_output)
    actual = parse_transition_output(transition_run.output_path)

    assert_metadata_matches(
        name=transition_run.name,
        expected=expected.metadata,
        actual=actual.metadata,
    )

    assert_metrics_match(
        name=transition_run.name,
        expected=expected.metrics,
        actual=actual.metrics,
    )

    assert actual.profile.keys() == expected.profile.keys(), (
        f"{transition_run.name} reported sequence positions changed"
    )

    for position, expected_row in expected.profile.items():
        actual_row = actual.profile[position]

        assert actual_row.base == expected_row.base, (
            f"{transition_run.name} base changed at position {position}: "
            f"expected {expected_row.base!r}, "
            f"actual {actual_row.base!r}"
        )

        assert_number_matches(
            name=f"{transition_run.name} position {position} P(x)",
            expected=expected_row.probability,
            actual=actual_row.probability,
        )

        assert (actual_row.energy is None) == (
            expected_row.energy is None
        ), (
            f"{transition_run.name} energy output changed at "
            f"position {position}"
        )

        if expected_row.energy is not None:
            assert actual_row.energy is not None

            assert_number_matches(
                name=f"{transition_run.name} position {position} G(x)",
                expected=expected_row.energy,
                actual=actual_row.energy,
            )
