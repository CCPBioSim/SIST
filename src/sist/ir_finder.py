"""Inverted repeat scoring for the cruciform and competition algorithms.

This module runs the external Inverted Repeat Finder (``irf``) tool
against a sequence, parses its HTML report, and scores each candidate
inverted repeat (IR) using the energy model in :mod:`sist.energetics`.
The result is the ``"start,length,energy,|"`` string consumed by
``qsidd -X``.

IRF download page: http://tandem.bu.edu/irf/irf.download.html
Reference: Zhabinskaya & Benham, Nucleic Acids Res, 41(21), 9610.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from Bio import SeqIO
from bs4 import BeautifulSoup

from sist._index import index_or_none
from sist.energetics import BASE_PAIR_COMPLEMENTS, CruciformEnergetics

_DNA_BASES = frozenset("ATCG")


class IRFinder:
    """Scores candidate inverted repeats for the cruciform/competition algorithms.

    Attributes:
        temperature: Temperature in Kelvin.
        shape: Either ``"linear"`` or ``"circular"``.
        energetics: Energy calculator for this temperature. Uses a fixed
            salt concentration independent of (and always overriding) any
            ``-i``/salt value passed to the ``qsidd`` binaries.
    """

    IRF_EXECUTABLE = "irf"
    IRF_MATCH = 2
    IRF_MISMATCH = 10
    IRF_DELTA = 10
    IRF_PM = 80
    IRF_PI = 10
    IRF_MIN_SCORE = 20
    IRF_MAX_LENGTH = 10000
    IRF_MAX_LOOP = 100
    IRF_MIN_LOOP = 3
    IRF_TIMEOUT_SECONDS = 120

    CRUCIFORM_SALT = 0.01

    def __init__(self, *, temperature: float, shape: str) -> None:
        """Initialise the finder.

        Args:
            temperature: Temperature in Kelvin.
            shape: Either ``"linear"`` or ``"circular"``.
        """
        self.temperature = temperature
        self.shape = shape
        self.energetics = CruciformEnergetics(
            temperature=temperature, salt=self.CRUCIFORM_SALT
        )

        self._seen: dict[tuple[int, int, int], int] = {}
        self._ir_string_parts: list[str] = ["1,0,10000,|"]

    def compute_cruciform_energy_string(self, file_path: str | Path) -> str:
        """Score every candidate inverted repeat in a sequence.

        Converts the input sequence to IRF's required format, runs IRF
        (twice, for circular sequences, to also score IRs spanning the
        origin), and returns the ``"start,length,energy,|..."`` string
        consumed by ``qsidd -X``.

        Args:
            file_path: Path to the input FASTA sequence file.

        Returns:
            The energy-scored inverted repeat string for the sequence.
        """
        one_line_path = self.convert_to_one_line(Path(file_path))
        sequence = self.read_sequence(one_line_path)

        cruciform_initiation_energy = self.energetics.cruciform_initiation_energy()

        self._score_reports(
            one_line_path,
            circular_pass=False,
            sequence=sequence,
            sequence_length=len(sequence),
            cruciform_initiation_energy=cruciform_initiation_energy,
        )

        if self.shape == "circular":
            sequence_length, shifted_sequence = self.sequence_length_and_shift(
                one_line_path, circular=True
            )

            circular_path = Path(f"circ.{one_line_path.name}")
            circular_path.write_text(
                f">{circular_path.name}\n{shifted_sequence}", encoding="utf-8"
            )

            self._score_reports(
                circular_path,
                circular_pass=True,
                sequence=shifted_sequence,
                sequence_length=sequence_length,
                cruciform_initiation_energy=cruciform_initiation_energy,
            )

        return "".join(self._ir_string_parts)

    @staticmethod
    def convert_to_one_line(sequence_path: Path) -> Path:
        """Convert a FASTA file into IRF's required single-line format.

        Uses Biopython to parse well-formed FASTA input. A file with no
        ``>`` header isn't valid FASTA, so Biopython rejects it; in that
        case, the first line is discarded entirely (not kept as sequence
        data) and a synthetic header is written instead.

        Args:
            sequence_path: Path to the input sequence file. The output
                file is written in the current working directory, named
                after this path's basename.

        Returns:
            Path to the newly written ``one_line.<name>`` file.
        """
        output_path = Path(f"one_line.{sequence_path.name}")

        try:
            records = list(SeqIO.parse(sequence_path, "fasta"))
        except ValueError:
            records = None

        if records:
            header = records[0].description
            sequence = str(records[0].seq)
        elif records == []:
            output_path.write_text("", encoding="utf-8")
            return output_path
        else:
            lines = sequence_path.read_text(encoding="utf-8").splitlines()
            header = output_path.name
            sequence = "".join(lines[1:])

        output_path.write_text(f">{header}\n{sequence}", encoding="utf-8")

        return output_path

    @staticmethod
    def read_sequence(one_line_path: Path) -> str:
        """Read the sequence out of a one-line FASTA file.

        Args:
            one_line_path: Path to a file in the single-line FASTA format
                produced by :meth:`convert_to_one_line`.

        Returns:
            The bare sequence, with no header line.
        """
        return str(SeqIO.read(one_line_path, "fasta").seq)

    @classmethod
    def sequence_length_and_shift(
        cls, one_line_path: Path, *, circular: bool
    ) -> tuple[int, str]:
        """Read a one-line sequence file and optionally shift it to the middle.

        Args:
            one_line_path: Path to a file in the single-line FASTA format
                produced by :meth:`convert_to_one_line`.
            circular: If True, also compute the sequence with its start
                position shifted to the middle (for circular re-indexing).

        Returns:
            A tuple of (sequence length, shifted sequence or "" if not
            circular).
        """
        sequence = cls.read_sequence(one_line_path)
        sequence_length = len(sequence)
        shifted_sequence = ""

        if circular:
            midpoint = sequence_length // 2
            shifted_sequence = sequence[midpoint:] + sequence[:midpoint]

        return sequence_length, shifted_sequence

    @staticmethod
    def parse_loop_header(line: str) -> tuple[int, int, int, int, int]:
        """Parse an IRF report ``Loop:`` line.

        Args:
            line: A line such as
                ``"    Indices: 2982--2992,2996--3006  Loop: 3  Score: 22"``.

        Returns:
            A tuple of (loop_length, loop_start, ir_start, ir_end, ir_length).
        """
        tokens = line.split()
        loop_length = int(tokens[3])

        left_range, right_range = tokens[1].split(",")
        ir_start, left_end = (int(value) for value in left_range.split("--"))
        _right_start, ir_end = (int(value) for value in right_range.split("--"))

        loop_start = left_end + 1
        ir_length = ir_end - ir_start + 1

        return loop_length, loop_start, ir_start, ir_end, ir_length

    @staticmethod
    def parse_arm(line: str) -> list[str]:
        """Parse one alignment token out of an IRF arm line.

        Args:
            line: A line such as ``"   2982 >> AAACCACCGCT >> 2992"`` or
                ``"   3006 << *********** << 2996"``.

        Returns:
            The arm's alignment tokens, one character each.
        """
        return list(line.split()[2])

    @staticmethod
    def _read_report_lines(report_path: Path) -> list[str]:
        """Read an IRF HTML report's text content, one entry per line.

        Args:
            report_path: Path to an IRF ``*.txt.html`` report.

        Returns:
            The report's text content (HTML tags stripped), split into
            lines.
        """
        soup = BeautifulSoup(report_path.read_text(encoding="utf-8"), "html.parser")
        pre = soup.find("pre")
        text = pre.get_text() if pre is not None else soup.get_text()

        return text.splitlines()

    @classmethod
    def _extend_loop_to_minimum(
        cls,
        loop_bases: list[str],
        left_arm: list[str],
        right_arm: list[str],
        *,
        loop_start: int,
    ) -> tuple[list[str], int, int]:
        """Borrow bases from the arms until the loop reaches ``IRF_MIN_LOOP``.

        Args:
            loop_bases: Bases in the loop, mutated in place.
            left_arm: Left-arm alignment tokens.
            right_arm: Right-arm alignment tokens, aligned with
                ``left_arm``.
            loop_start: Loop start position, before extension.

        Returns:
            A tuple of (extended loop bases, updated loop_start,
            arm bases consumed by the extension).
        """
        shorten_arm = 0
        index = -1

        while len(loop_bases) < cls.IRF_MIN_LOOP:
            left_value = index_or_none(left_arm, index)

            if left_value is not None and left_value in _DNA_BASES:
                right_value = index_or_none(right_arm, index)

                if right_value == "-":
                    loop_bases.append(left_value)
                elif right_value == "*":
                    loop_bases.append(left_value)
                    loop_bases.append(BASE_PAIR_COMPLEMENTS[left_value])
                else:
                    assert right_value is not None
                    loop_bases.append(left_value)
                    loop_bases.append(right_value)
            else:
                loop_bases.append(index_or_none(right_arm, index))  # type: ignore[arg-type]
                shorten_arm += 1

            loop_start -= 1
            shorten_arm += 1
            index -= 1

        return loop_bases, loop_start, shorten_arm

    def _score_reports(
        self,
        sequence_file: Path,
        *,
        circular_pass: bool,
        sequence: str,
        sequence_length: int,
        cruciform_initiation_energy: float,
    ) -> None:
        """Run IRF against ``sequence_file`` and append scored IRs to the result.

        Args:
            sequence_file: One-line FASTA file to run IRF against.
            circular_pass: True when scoring the circularly shifted
                sequence.
            sequence: The bare sequence read from ``sequence_file``.
            sequence_length: Length of the (unshifted) original sequence.
            cruciform_initiation_energy: Precomputed Ecr for this
                temperature.
        """
        subprocess.run(
            [
                self.IRF_EXECUTABLE,
                str(sequence_file),
                str(self.IRF_MATCH),
                str(self.IRF_MISMATCH),
                str(self.IRF_DELTA),
                str(self.IRF_PM),
                str(self.IRF_PI),
                str(self.IRF_MIN_SCORE),
                str(self.IRF_MAX_LENGTH),
                str(self.IRF_MAX_LOOP),
            ],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=self.IRF_TIMEOUT_SECONDS,
        )

        report_prefix = (
            f"{sequence_file}.{self.IRF_MATCH}.{self.IRF_MISMATCH}.{self.IRF_DELTA}."
            f"{self.IRF_PM}.{self.IRF_PI}.{self.IRF_MIN_SCORE}.{self.IRF_MAX_LENGTH}."
            f"{self.IRF_MAX_LOOP}"
        )

        report_index = 1

        while (
            report_path := Path(f"{report_prefix}.{report_index}.txt.html")
        ).is_file():
            report_index += 1

            reading = False
            skip = False
            left_arm: list[str] = []
            right_arm: list[str] = []
            loop_length = loop_start = ir_start = ir_end = ir_length = 0

            for line in self._read_report_lines(report_path):
                if "frequent" in line:
                    reading = False

                if "Done" in line:
                    break

                if "Loop:" in line:
                    skip = False
                    loop_length, loop_start, ir_start, ir_end, ir_length = (
                        self.parse_loop_header(line)
                    )

                    if loop_length > self.IRF_MAX_LOOP:
                        continue

                    self._seen[(loop_start, loop_length, ir_start)] = ir_length
                    reading = True

                if not reading:
                    continue

                if ">>" in line and "LF" not in line:
                    left_arm.extend(self.parse_arm(line))

                if "<<" in line and "RF" not in line:
                    right_arm.extend(self.parse_arm(line))

                if "Statistics" in line:
                    loop_bases = list(sequence[loop_start : loop_start + loop_length])

                    loop_bases, loop_start, shorten_arm = self._extend_loop_to_minimum(
                        loop_bases, left_arm, right_arm, loop_start=loop_start
                    )
                    loop_length = len(loop_bases)

                    if circular_pass:
                        half = sequence_length // 2

                        if loop_start > half:
                            offset = half + (sequence_length % 2)
                            loop_start -= offset
                            ir_start -= offset
                        else:
                            loop_start += half
                            ir_start += half

                        seen_length = self._seen.get(
                            (loop_start, loop_length, ir_start)
                        )

                        if seen_length is not None and seen_length == ir_length:
                            skip = True
                            left_arm.clear()
                            right_arm.clear()

                    if skip:
                        continue

                    loop_free_energy = self.energetics.loop_energy(loop_bases)

                    arm_length = len(left_arm) - shorten_arm
                    left_arm = left_arm[:arm_length]
                    right_arm = right_arm[:arm_length]

                    imperfection_energies = self.energetics.imperfection_energies(
                        right_arm, left_arm
                    )
                    imperfection_energies.reverse()

                    total_energy = cruciform_initiation_energy + loop_free_energy
                    ir_position = loop_start
                    ir_length_running = loop_length
                    extension_count = 0

                    for position in range(arm_length):
                        total_energy += 2 * imperfection_energies[position]

                        if left_arm[arm_length - 1 - position] != "-":
                            extension_count += 1
                            ir_position = loop_start - extension_count
                            ir_length_running = loop_length + 2 * extension_count

                        self._ir_string_parts.append(
                            f"{ir_position},{ir_length_running},{total_energy},|"
                        )

                    left_arm = []
                    right_arm = []
