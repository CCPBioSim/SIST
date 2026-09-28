"""Unit tests for sist.energetics, hand-verified against the reference paper.

Expected values are computed independently from the published formulas
(not by re-calling the module under test), so these tests catch
transcription errors.
"""

from __future__ import annotations

import math

import pytest

from sist import energetics

TEMPERATURE = 310.0
SALT = 0.01


@pytest.fixture
def calculator() -> energetics.CruciformEnergetics:
    return energetics.CruciformEnergetics(temperature=TEMPERATURE, salt=SALT)


def _reference_melting_energy(base: str, *, temperature: float, salt: float) -> float:
    at_melting_temperature = 354.65 + 16.6 * math.log10(salt)
    gc_melting_temperature = at_melting_temperature + 41
    at_energy = 7.2464 * (1 - temperature / at_melting_temperature)
    gc_energy = 9.0172 * (1 - temperature / gc_melting_temperature)
    return at_energy if base in ("a", "t", "A", "T") else gc_energy


@pytest.mark.parametrize("base", ["A", "a", "T", "t"])
def test_melting_energy_at_bases(
    calculator: energetics.CruciformEnergetics, base: str
) -> None:
    expected = _reference_melting_energy(base, temperature=TEMPERATURE, salt=SALT)

    actual = calculator.melting_energy(base)

    assert actual == pytest.approx(expected)


@pytest.mark.parametrize("base", ["G", "g", "C", "c"])
def test_melting_energy_gc_bases(
    calculator: energetics.CruciformEnergetics, base: str
) -> None:
    expected = _reference_melting_energy(base, temperature=TEMPERATURE, salt=SALT)

    actual = calculator.melting_energy(base)

    assert actual == pytest.approx(expected)


def test_gas_constant_rt(calculator: energetics.CruciformEnergetics) -> None:
    assert calculator.gas_constant_rt == pytest.approx(1.9872 * TEMPERATURE / 1000.0)


def test_loop_energy_matches_reference_formula(
    calculator: energetics.CruciformEnergetics,
) -> None:
    loop = ["A", "T", "G", "C", "A"]

    expected = sum(
        _reference_melting_energy(base, temperature=TEMPERATURE, salt=SALT)
        for base in loop
    )
    expected += 2 * 2.44 * (1.9872 * TEMPERATURE / 1000.0) * math.log(len(loop))

    actual = calculator.loop_energy(loop)

    assert actual == pytest.approx(expected)


def test_mismatch_energy_matches_reference_tables() -> None:
    triplet1, triplet2 = "GAA", "CAT"

    e_wc = (
        energetics.WATSON_CRICK_ENERGIES["GA"] + energetics.WATSON_CRICK_ENERGIES["AA"]
    )
    e_miss = (
        energetics.MISMATCH_ENERGIES["GA.CA"] + energetics.MISMATCH_ENERGIES["TA.AA"]
    )
    expected = e_miss - e_wc

    actual = energetics.CruciformEnergetics.mismatch_energy(triplet1, triplet2)

    assert actual == pytest.approx(expected)


def test_cruciform_initiation_energy_matches_reference_formula(
    calculator: energetics.CruciformEnergetics,
) -> None:
    rt = 1.9872 * TEMPERATURE / 1000.0
    expected = (
        192.5
        - TEMPERATURE * 0.565
        - 4 * _reference_melting_energy("A", temperature=TEMPERATURE, salt=SALT)
        - 2 * 2.44 * rt * math.log(4)
    )

    actual = calculator.cruciform_initiation_energy()

    assert actual == pytest.approx(expected)


def test_imperfection_energies_requires_star_at_start(
    calculator: energetics.CruciformEnergetics,
) -> None:
    with pytest.raises(ValueError, match="star"):
        calculator.imperfection_energies(["A", "*"], ["A", "A"])


def test_imperfection_energies_perfect_match_is_all_zero(
    calculator: energetics.CruciformEnergetics,
) -> None:
    right_arm = ["*", "*", "*", "*"]
    left_arm = ["A", "T", "G", "C"]

    energies = calculator.imperfection_energies(right_arm, left_arm)

    assert energies == [0.0, 0.0, 0.0, 0.0]


def test_imperfection_energies_bulge_on_left_arm(
    calculator: energetics.CruciformEnergetics,
) -> None:
    # Right arm has a gap at position 1 -> bulge on the left arm base "T".
    right_arm = ["*", "-", "*"]
    left_arm = ["A", "T", "G"]

    energies = calculator.imperfection_energies(right_arm, left_arm)

    expected_bulge = energetics.IMPERFECTION_ENERGY + _reference_melting_energy(
        "T", temperature=TEMPERATURE, salt=SALT
    )

    assert energies[1] == pytest.approx(expected_bulge)
    assert energies[0] == pytest.approx(0.0)
    assert energies[2] == pytest.approx(0.0)


def test_imperfection_energies_bulge_on_right_arm(
    calculator: energetics.CruciformEnergetics,
) -> None:
    # Left arm has a gap at position 1 -> bulge on the right arm base "C".
    right_arm = ["*", "C", "*"]
    left_arm = ["A", "-", "G"]

    energies = calculator.imperfection_energies(right_arm, left_arm)

    expected_bulge = energetics.IMPERFECTION_ENERGY + _reference_melting_energy(
        "C", temperature=TEMPERATURE, salt=SALT
    )

    assert energies[1] == pytest.approx(expected_bulge)


def test_imperfection_energies_internal_loop(
    calculator: energetics.CruciformEnergetics,
) -> None:
    # Two adjacent open positions (right arm not '*' on both sides of index 1).
    right_arm = ["*", "-", "A", "-", "*"]
    left_arm = ["A", "T", "T", "T", "G"]

    energies = calculator.imperfection_energies(right_arm, left_arm)

    expected_internal_loop = (
        2 * energetics.IMPERFECTION_ENERGY
        + _reference_melting_energy("A", temperature=TEMPERATURE, salt=SALT)
        + _reference_melting_energy("T", temperature=TEMPERATURE, salt=SALT)
    )

    assert energies[2] == pytest.approx(expected_internal_loop)


def test_imperfection_energies_mismatch(
    calculator: energetics.CruciformEnergetics,
) -> None:
    # A single mismatch flanked by perfect matches on both sides.
    right_arm = ["*", "A", "*"]
    left_arm = ["G", "G", "T"]

    energies = calculator.imperfection_energies(right_arm, left_arm)

    expected = energetics.CruciformEnergetics.mismatch_energy("GGT", "CAA")

    assert energies[1] == pytest.approx(expected)
