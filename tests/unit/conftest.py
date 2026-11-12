"""Shared fixtures for the unit tests."""

from __future__ import annotations

import math
from typing import Any

import pytest

SEQUENCE = "TTCTC"
SEQUENCE_SHA256 = "2f02431f1d80b9dd4c4da090b4aae94406b6129aa7f78d95f6e4f19b3a30c3e2"


def _base_result(algorithm: str, transitions: list[str]) -> dict[str, Any]:
    """Build the parts every SIST result shares, for a 5 bp sequence.

    Args:
        algorithm: Algorithm flag, such as ``"M"``.
        transitions: Transitions the algorithm calculates.

    Returns:
        A result with everything except the summary and profile filled in.
    """
    return {
        "schema_version": "1.0.0",
        "status": "success",
        "provenance": {
            "software": {
                "name": "SIST",
                "version": "1.1.0",
                "git_sha": None,
                "engine": "trans_compete" if algorithm == "A" else "trans_three",
                "irf_version": "3.09" if algorithm in "CA" else None,
            },
            "environment": {"python": "3.13.13", "platform": "Linux-x86_64"},
            "command": ["sist", "-a", algorithm, "-f", "sequence.fa"],
            "started_at": "2026-11-06T10:00:00Z",
            "finished_at": "2026-11-06T10:00:02Z",
            "wall_time_s": 2.1,
        },
        "calculation": {"algorithm": algorithm, "transitions": transitions},
        "input": {
            "sequence_file": "sequence.fa",
            "sequence_name": None,
            "sequence_length_bp": len(SEQUENCE),
            "sequence_sha256": SEQUENCE_SHA256,
            "molecule_type": "linear",
        },
        "parameters": {
            "temperature_k": 310.0,
            "superhelical_density": -0.06,
            "salt_concentration_m": 0.01,
            "energy_threshold_kcal_per_mol": 12.0,
            "melting_energetics": "copolymeric" if "melting" in transitions else None,
        },
        "derived_parameters": {
            "linking_difference_turns": -25.448076923076918,
            "torsional_stiffness_kcal_bp_per_rad2": 3.6,
            "nucleation_energy_kcal_per_mol": {
                transition: 10.84 for transition in transitions
            },
            "bat_kcal_per_mol": 0.2581156634002174,
            "bgc_kcal_per_mol": 1.3048755414539936,
        },
        "warnings": [],
        "error": None,
    }


def _summary(transitions: list[str], partition_function: float) -> dict[str, Any]:
    """Build a summary with the same averages for every transition.

    Args:
        transitions: Transitions to summarise.
        partition_function: Partition function to record.

    Returns:
        A summary whose state counts add up and whose log matches.
    """
    return {
        "ln_partition_function": math.log(partition_function),
        "partition_function": partition_function,
        "average_number_of_runs": 1.94048,
        "state_counts": {
            "by_number_of_runs": [10652, 4908539, 2873551, 0],
            "total": 7792742,
        },
        "transitions": {
            transition: {
                "probability": 0.7,
                "expected_runs": 1.02,
                "expected_transformed_bp": 24.7,
            }
            for transition in transitions
        },
    }


@pytest.fixture
def melting_result() -> dict[str, Any]:
    """A valid melting (M) result for a 5 bp sequence."""
    result = _base_result("M", ["melting"])

    result["summary"] = _summary(["melting"], 2.22693e-60)
    result["profile"] = {
        "base": SEQUENCE,
        "probability": {"melting": [2.4e-07, 3.3e-07, 0.0, 5.2e-07, 5.9e-07]},
        "free_energy_kcal_per_mol": {"melting": [11.06, 11.08, None, 11.03, 11.04]},
    }

    return result


@pytest.fixture
def competition_result() -> dict[str, Any]:
    """A valid competition (A) result for a 5 bp sequence."""
    transitions = ["melting", "z_dna", "cruciform"]
    result = _base_result("A", transitions)
    result["summary"] = _summary(transitions, 1.4748e-54)
    result["summary"]["pairwise"] = {
        name: {
            first: {second: 0.69 for second in transitions if second != first}
            for first in transitions
        }
        for name in (
            "joint_probability",
            "conditional_probability",
            "conditional_probability_given_absent",
        )
    }
    result["profile"] = {
        "base": SEQUENCE,
        "probability": {
            transition: [6.6e-08, 8.1e-08, 6.4e-08, 1.0e-07, 9.8e-08]
            for transition in transitions
        },
    }
    result["inverted_repeats"] = [
        {"start": 2, "length_bp": 3, "energy_kcal_per_mol": 17.27}
    ]
    result["parameters"]["cruciform_salt_concentration_m"] = 0.01

    return result


@pytest.fixture
def failed_result() -> dict[str, Any]:
    """A valid failed competition (A) run, where IRF was not installed."""
    result = _base_result("A", ["melting", "z_dna", "cruciform"])

    result["status"] = "failed"
    result["derived_parameters"] = None
    result["summary"] = None
    result["profile"] = None
    result["error"] = {
        "code": "irf_not_found",
        "message": "The irf executable was not found on PATH.",
        "stage": "ir_finder",
    }

    return result
