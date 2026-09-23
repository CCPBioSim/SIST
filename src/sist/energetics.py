"""Thermodynamic energy calculations for cruciform extrusion candidates.

Used to score candidate inverted repeats (IRs) reported by the Inverted
Repeat Finder (IRF) tool before they are handed to the ``qsidd``
cruciform/competition binaries.

Reference: Zhabinskaya & Benham, Nucleic Acids Res, 41(21), 9610.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from sist._index import index_or_none

BASE_PAIR_COMPLEMENTS: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}

WATSON_CRICK_ENERGIES: dict[str, float] = {
    "AA": -1.0,
    "TT": -1.0,
    "AT": -0.88,
    "TA": -0.58,
    "CA": -1.45,
    "AC": -1.45,
    "GT": -1.44,
    "TG": -1.44,
    "CT": -1.28,
    "TC": -1.28,
    "GA": -1.3,
    "AG": -1.3,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.84,
    "CC": -1.84,
}

MISMATCH_ENERGIES: dict[str, float] = {
    "GA.CA": 0.17,
    "GA.CC": 0.81,
    "GA.CG": -0.25,
    "GC.CA": 0.47,
    "GC.CC": 0.79,
    "GC.CT": 0.62,
    "GG.CA": -0.52,
    "GG.CG": -1.11,
    "GG.CT": 0.08,
    "GT.CC": 0.98,
    "GT.CG": -0.59,
    "GT.CT": 0.45,
    "CA.GA": 0.43,
    "CA.GC": 0.75,
    "CA.GG": 0.03,
    "CC.GA": 0.79,
    "CC.GC": 0.70,
    "CC.GT": 0.62,
    "CG.GA": 0.11,
    "CG.GG": -0.11,
    "CG.GT": -0.47,
    "CT.GC": 0.40,
    "CT.GG": -0.32,
    "CT.GT": -0.12,
    "AA.TA": 0.61,
    "AA.TC": 0.88,
    "AA.TG": 0.14,
    "AC.TA": 0.77,
    "AC.TC": 1.33,
    "AC.TT": 0.64,
    "AG.TA": 0.02,
    "AG.TG": -0.13,
    "AG.TT": 0.71,
    "AT.TC": 0.73,
    "AT.TG": 0.07,
    "AT.TT": 0.69,
    "TA.AA": 0.69,
    "TA.AC": 0.92,
    "TA.AG": 0.42,
    "TC.AA": 1.33,
    "TC.AC": 1.05,
    "TC.AT": 0.97,
    "TG.AA": 0.74,
    "TG.AG": 0.44,
    "TG.AT": 0.43,
    "TT.AC": 0.75,
    "TT.AG": 0.34,
    "TT.AT": 0.68,
}

IMPERFECTION_ENERGY = 2.42


class CruciformEnergetics:
    """Energy calculations for cruciform extrusion candidates.

    All calculations are evaluated at a fixed temperature and salt
    concentration, set once when the calculator is constructed.

    Attributes:
        temperature: Temperature in Kelvin.
        salt: Salt (ionic strength) concentration in mol/L.
    """

    def __init__(self, *, temperature: float, salt: float) -> None:
        """Initialise the calculator.

        Args:
            temperature: Temperature in Kelvin.
            salt: Salt (ionic strength) concentration in mol/L.
        """
        self.temperature = temperature
        self.salt = salt

    @property
    def gas_constant_rt(self) -> float:
        """RT in kcal/mol at this calculator's temperature."""
        return 1.9872 * self.temperature / 1000.0

    def melting_energy(self, base: str) -> float:
        """Compute the melting energy of a single base pair.

        Args:
            base: The base at this position. ``A``/``T`` (case-insensitive)
                use the AT melting energy; anything else uses the GC
                melting energy.

        Returns:
            Melting energy in kcal/mol.
        """
        at_melting_temperature = 354.65 + 16.6 * math.log10(self.salt)
        gc_melting_temperature = at_melting_temperature + 41

        at_energy = 7.2464 * (1 - self.temperature / at_melting_temperature)
        gc_energy = 9.0172 * (1 - self.temperature / gc_melting_temperature)

        if base in ("a", "t", "A", "T"):
            return at_energy

        return gc_energy

    def loop_energy(self, loop_sequence: Sequence[str]) -> float:
        """Compute the free energy of a cruciform loop.

        Args:
            loop_sequence: Bases making up the loop.

        Returns:
            Loop free energy in kcal/mol, including the loop entropy term.
        """
        energy = sum(self.melting_energy(base) for base in loop_sequence)
        energy += 2 * 2.44 * self.gas_constant_rt * math.log(len(loop_sequence))

        return energy

    @staticmethod
    def mismatch_energy(triplet1: str, triplet2: str) -> float:
        """Compute the energy penalty of a base-pair mismatch.

        Args:
            triplet1: Three consecutive bases on one strand, centred on the
                mismatch.
            triplet2: The corresponding three bases on the opposite strand.

        Returns:
            Mismatch energy penalty in kcal/mol.
        """
        watson_crick_energy = (
            WATSON_CRICK_ENERGIES[triplet1[0] + triplet1[1]]
            + WATSON_CRICK_ENERGIES[triplet1[1] + triplet1[2]]
        )

        stack1 = f"{triplet1[0]}{triplet1[1]}.{triplet2[0]}{triplet2[1]}"
        stack2 = f"{triplet2[2]}{triplet2[1]}.{triplet1[2]}{triplet1[1]}"
        mismatch_stack_energy = MISMATCH_ENERGIES[stack1] + MISMATCH_ENERGIES[stack2]

        return mismatch_stack_energy - watson_crick_energy

    def imperfection_energies(
        self, right_arm: Sequence[str], left_arm: Sequence[str]
    ) -> list[float]:
        """Compute per-position imperfection energies along an IR arm alignment.

        ``right_arm`` and ``left_arm`` are aligned, equal-length sequences of
        single-character tokens describing one arm of a candidate inverted
        repeat: a base letter for a bulge, ``"-"`` for a gap, or (on
        ``right_arm``) ``"*"`` for a perfectly matched position.

        Args:
            right_arm: Right-arm alignment tokens, starting with ``"*"``.
            left_arm: Left-arm alignment tokens, aligned with ``right_arm``.

        Returns:
            One imperfection energy (kcal/mol) per position, ``0.0`` where
            there is no imperfection.

        Raises:
            ValueError: If ``right_arm`` does not start with ``"*"``.
        """
        if right_arm[0] != "*":
            raise ValueError("Doesn't start with a star")

        energies = [0.0] * len(left_arm)

        for index, left_base in enumerate(left_arm):
            if left_base == "-":
                energies[index] = IMPERFECTION_ENERGY + self.melting_energy(
                    right_arm[index]
                )

        for index, right_base in enumerate(right_arm):
            if right_base == "-":
                energies[index] = IMPERFECTION_ENERGY + self.melting_energy(
                    left_arm[index]
                )
                continue

            if right_base in ("*", "-"):
                continue

            if left_arm[index] == "-":
                continue

            previous_right = index_or_none(right_arm, index - 1)
            next_right = index_or_none(right_arm, index + 1)

            if previous_right != "*" or next_right != "*":
                energies[index] = (
                    2 * IMPERFECTION_ENERGY
                    + self.melting_energy(right_base)
                    + self.melting_energy(left_arm[index])
                )
                continue

            previous_left = index_or_none(left_arm, index - 1)
            next_left = index_or_none(left_arm, index + 1)
            assert previous_left is not None
            assert next_left is not None

            if previous_left != "-" and next_left != "-":
                triplet1 = f"{previous_left}{left_arm[index]}{next_left}"
                triplet2 = (
                    f"{BASE_PAIR_COMPLEMENTS[previous_left]}"
                    f"{right_base}"
                    f"{BASE_PAIR_COMPLEMENTS[next_left]}"
                )
                energies[index] = self.mismatch_energy(triplet1, triplet2)

        return energies

    def cruciform_initiation_energy(self) -> float:
        """Compute the cruciform initiation energy, Ecr.

        Returns:
            Cruciform initiation energy in kcal/mol.
        """
        return (
            192.5
            - self.temperature * 0.565
            - 4 * self.melting_energy("A")
            - 2 * 2.44 * self.gas_constant_rt * math.log(4)
        )
