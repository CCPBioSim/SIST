"""Orchestration that dispatches a SIST calculation to the qsidd binaries.

Melting (M) and Z-DNA (Z) calculations run the single-transition
``qsidd`` binary directly; cruciform (C) and competition (A) calculations
first score candidate inverted repeats via :mod:`sist.ir_finder`, then run
the single-transition or competition ``qsidd`` binary respectively.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from sist.binaries import BinaryLocator
from sist.ir_finder import IRFinder


class SistRunner:
    """Coordinate a single SIST calculation: build, run, and report on it.

    Attributes:
        file: Path to the input FASTA sequence file.
        algorithm: One of ``"M"`` (melting), ``"Z"`` (Z-DNA),
            ``"C"`` (cruciform), or ``"A"`` (competition).
        temperature: Temperature in Kelvin.
        superhelical_density: Superhelical density (sigma).
        salt: Salt (ionic strength) concentration in mol/L.
        energy_threshold: Energy threshold for the transition search.
        circular: Whether the input molecule is circular (default linear).
        nearest_neighbor: Use nearest-neighbor melting energetics.
            Ignored for the Z-DNA and cruciform algorithms.
        print_base_pair: Print the base pair for each position.
        print_parameters: Print algorithm parameters.
        print_ensemble_average: Print ensemble average results.
        output_file: Path to write output to, or None to print to stdout.
    """

    def __init__(
        self,
        *,
        file: str,
        algorithm: str,
        temperature: float,
        superhelical_density: float,
        salt: float,
        energy_threshold: float,
        circular: bool,
        nearest_neighbor: bool,
        print_base_pair: bool,
        print_parameters: bool,
        print_ensemble_average: bool,
        output_file: str | None,
    ) -> None:
        """Initialise the runner with a single calculation's configuration."""
        self.file = file
        self.algorithm = algorithm
        self.temperature = temperature
        self.superhelical_density = superhelical_density
        self.salt = salt
        self.energy_threshold = energy_threshold
        self.circular = circular
        self.nearest_neighbor = nearest_neighbor
        self.print_base_pair = print_base_pair
        self.print_parameters = print_parameters
        self.print_ensemble_average = print_ensemble_average
        self.output_file = output_file

        self._binaries = BinaryLocator()

    @property
    def shape(self) -> str:
        """Molecule shape as the qsidd binaries expect it: linear or circular."""
        return "circular" if self.circular else "linear"

    def run(self) -> None:
        """Run the calculation and write or print its output.

        Raises:
            ValueError: If ``algorithm`` is not one of M, Z, C, or A.
            subprocess.CalledProcessError: If the qsidd binary exits
                non-zero.
        """
        command = self._build_command()
        result = subprocess.run(command, capture_output=True, text=True, check=True)

        if self.output_file:
            Path(self.output_file).write_text(result.stdout, encoding="utf-8")
        else:
            print(result.stdout, end="")

    def _build_command(self) -> list[str]:
        """Build the qsidd command line for this calculation's algorithm."""
        shared_parameters = self._shared_parameters()

        if self.algorithm == "M":
            return [
                str(self._binaries.trans_three_binary()),
                *self._flags(nearest_neighbor=self.nearest_neighbor),
                *shared_parameters,
                "-f",
                self.file,
            ]

        if self.algorithm == "Z":
            return [
                str(self._binaries.trans_three_binary()),
                *self._flags(nearest_neighbor=False),
                *shared_parameters,
                "-Z",
                "-f",
                self.file,
            ]

        if self.algorithm == "C":
            ir_string = self._cruciform_energy_string()

            return [
                str(self._binaries.trans_three_binary()),
                *self._flags(nearest_neighbor=False),
                *shared_parameters,
                "-C",
                "-X",
                ir_string,
                "-f",
                self.file,
            ]

        if self.algorithm == "A":
            ir_string = self._cruciform_energy_string()

            return [
                str(self._binaries.trans_compete_binary()),
                *self._flags(nearest_neighbor=self.nearest_neighbor),
                *shared_parameters,
                "-X",
                ir_string,
                "-f",
                self.file,
            ]

        raise ValueError(
            f"Unknown algorithm: {self.algorithm!r}. Expected one of M, Z, C, A."
        )

    def _cruciform_energy_string(self) -> str:
        """Score candidate inverted repeats for the C/A algorithms."""
        finder = IRFinder(temperature=self.temperature, shape=self.shape)

        return finder.compute_cruciform_energy_string(self.file)

    def _flags(self, *, nearest_neighbor: bool) -> list[str]:
        """Build the qsidd boolean flags."""
        flags = []

        if self.print_base_pair:
            flags.append("-b")
        if self.print_parameters:
            flags.append("-p")
        if self.print_ensemble_average:
            flags.append("-r")
        if nearest_neighbor:
            flags.append("-n")
        if self.circular:
            flags.append("-c")

        return flags

    def _shared_parameters(self) -> list[str]:
        """Build the qsidd numeric parameter flags shared by every algorithm."""
        return [
            "-T",
            str(self.temperature),
            "-s",
            str(self.superhelical_density),
            "-i",
            str(self.salt),
            "-t",
            str(self.energy_threshold),
        ]
