"""Unit tests for sist.runner, using stand-ins for IRF and qsidd.

These cover the runner's own jobs: building qsidd commands, choosing between
text and JSON output, and turning failures into failed results. What goes into
a result is covered by the sist.results tests.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from sist.results.diagnostics import ResultError
from sist.runner import ScoredRepeats, SistRunner, SubprocessQsidd

QSIDD_PATH = Path(__file__).resolve().parent / "data" / "qsidd"


class FakeQsidd:
    """Stand-in for qsidd that prints a saved output and records each command.

    Attributes:
        output_name: File name of the saved qsidd output to print.
        returncode: Exit code to report.
        stderr: Standard error to report.
        replace: Text to swap in the saved output, as ``(old, new)``.
        commands: Every command line it was asked to run.
    """

    def __init__(
        self,
        output_name: str,
        *,
        returncode: int = 0,
        stderr: str = "",
        replace: tuple[str, str] = ("", ""),
    ) -> None:
        """Set up the stand-in.

        Args:
            output_name: File name of the saved qsidd output to print.
            returncode: Exit code to report.
            stderr: Standard error to report.
            replace: Text to swap in the saved output, as ``(old, new)``.
        """
        self.output_name = output_name
        self.returncode = returncode
        self.stderr = stderr
        self.replace = replace
        self.commands: list[list[str]] = []

    def run(self, command: list[str]) -> subprocess.CompletedProcess[str]:
        """Record the command and report the saved output.

        Args:
            command: The command line to run.

        Returns:
            The saved output as a finished process.
        """
        self.commands.append(command)
        stdout = (QSIDD_PATH / self.output_name).read_text(encoding="utf-8")

        if self.replace[0]:
            stdout = stdout.replace(*self.replace)

        return subprocess.CompletedProcess(
            command, self.returncode, stdout, self.stderr
        )


class FakeScorer:
    """Stand-in for the IR finder.

    Attributes:
        missing_irf: Whether to behave as if IRF is not installed.
    """

    def __init__(self, *, missing_irf: bool = False) -> None:
        """Set up the stand-in.

        Args:
            missing_irf: Whether to behave as if IRF is not installed.
        """
        self.missing_irf = missing_irf

    def score(self, file: str) -> ScoredRepeats:
        """Return the scored inverted repeats saved for the fragment.

        Args:
            file: Path to the sequence file.

        Returns:
            The saved inverted repeats.

        Raises:
            FileNotFoundError: If behaving as if IRF is not installed.
        """
        if self.missing_irf:
            raise FileNotFoundError("irf")

        ir_string = (QSIDD_PATH / "fragment.ir").read_text(encoding="utf-8").strip()

        return ScoredRepeats(ir_string, "3.09")


@pytest.fixture
def directory(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point SIST at placeholder qsidd binaries."""
    placeholder = tmp_path / "qsidd"
    placeholder.write_text("", encoding="utf-8")
    monkeypatch.setenv("SIST_TRANS_THREE_BIN", str(placeholder))
    monkeypatch.setenv("SIST_TRANS_COMPETE_BIN", str(placeholder))

    return tmp_path


def _runner(directory: Path, qsidd: FakeQsidd, **overrides: Any) -> SistRunner:
    """Make a runner for the 20 bp fragment, using the stand-ins.

    Args:
        directory: Directory holding the sequence and the output files.
        qsidd: The stand-in for qsidd.
        **overrides: Settings to change from the defaults.

    Returns:
        The runner.
    """
    sequence_path = directory / "fragment.fa"
    shutil.copy2(QSIDD_PATH / "fragment.fa", sequence_path)

    settings: dict[str, Any] = {
        "file": str(sequence_path),
        "algorithm": "M",
        "temperature": 310.0,
        "superhelical_density": 0.06,
        "salt": 0.01,
        "energy_threshold": 12.0,
        "circular": False,
        "nearest_neighbor": False,
        "print_base_pair": False,
        "print_parameters": False,
        "print_ensemble_average": False,
        "output_file": None,
        "json_file": str(directory / "result.json"),
        "command": ["sist", "-a", "M", "-f", "fragment.fa"],
        "ir_scorer": FakeScorer(),
        "qsidd": qsidd,
    }
    settings.update(overrides)

    return SistRunner(**settings)


def _read_result(directory: Path) -> dict[str, Any]:
    """Read the result the runner wrote.

    Args:
        directory: Directory the result was written to.

    Returns:
        The result.
    """
    return json.loads((directory / "result.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize(
    ("flag", "output_name", "algorithm_flags"),
    [
        ("M", "melting.txt", []),
        ("Z", "z-dna.txt", ["-Z"]),
        ("C", "cruciform.txt", ["-C", "-X"]),
        ("A", "competition.txt", ["-X"]),
    ],
)
def test_run_builds_the_command_from_the_algorithm(
    directory: Path, flag: str, output_name: str, algorithm_flags: list[str]
) -> None:
    qsidd = FakeQsidd(output_name)

    _runner(directory, qsidd, algorithm=flag).run()

    command = qsidd.commands[0]

    assert command[1:5] == ["-b", "-p", "-r", "-P"]
    assert [part for part in command if part in ("-Z", "-C", "-X")] == algorithm_flags
    assert _read_result(directory)["calculation"]["algorithm"] == flag


def test_run_passes_nearest_neighbor_only_where_supported(directory: Path) -> None:
    qsidd = FakeQsidd("z-dna.txt")

    _runner(directory, qsidd, algorithm="Z", nearest_neighbor=True).run()

    assert "-n" not in qsidd.commands[0]


def test_run_records_the_command_line(directory: Path) -> None:
    _runner(directory, FakeQsidd("melting.txt")).run()

    command = _read_result(directory)["provenance"]["command"]

    assert command == ["sist", "-a", "M", "-f", "fragment.fa"]


def test_run_with_only_json_prints_no_text(
    directory: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    _runner(directory, FakeQsidd("melting.txt")).run()

    assert capsys.readouterr().out == ""


def test_run_writes_text_from_a_second_normal_run(directory: Path) -> None:
    qsidd = FakeQsidd("melting.txt")
    output_path = directory / "result.txt"

    _runner(directory, qsidd, output_file=str(output_path), print_base_pair=True).run()

    assert output_path.read_text(encoding="utf-8").startswith("DNA sequence file")
    assert "-P" in qsidd.commands[0]
    assert "-P" not in qsidd.commands[1]
    assert "-b" in qsidd.commands[1]


def test_run_without_json_prints_the_text_output(
    directory: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    _runner(directory, FakeQsidd("melting.txt"), json_file=None).run()

    expected = (QSIDD_PATH / "melting.txt").read_text(encoding="utf-8")

    assert capsys.readouterr().out == expected
    assert not (directory / "result.json").exists()


def test_run_without_json_raises_when_qsidd_fails(directory: Path) -> None:
    runner = _runner(directory, FakeQsidd("melting.txt", returncode=3), json_file=None)

    with pytest.raises(subprocess.CalledProcessError):
        runner.run()


def test_run_rejects_an_unknown_algorithm(directory: Path) -> None:
    with pytest.raises(ValueError):
        _runner(directory, FakeQsidd("melting.txt"), algorithm="X").run()


@pytest.mark.parametrize(
    ("overrides", "qsidd", "code"),
    [
        ({}, FakeQsidd("melting.txt", returncode=3, stderr="broken"), "qsidd_failed"),
        ({}, FakeQsidd("fragment.fa"), "internal_error"),
        (
            {},
            FakeQsidd(
                "melting.txt",
                replace=("Average number of runs = ", "Average number of runs = -nan "),
            ),
            "internal_error",
        ),
        (
            {"algorithm": "C", "ir_scorer": FakeScorer(missing_irf=True)},
            FakeQsidd("cruciform.txt"),
            "irf_not_found",
        ),
    ],
    ids=["qsidd-fails", "unreadable-output", "nan-in-output", "irf-missing"],
)
def test_run_writes_a_failed_result(
    directory: Path, overrides: dict[str, Any], qsidd: FakeQsidd, code: str
) -> None:
    with pytest.raises(ResultError) as raised:
        _runner(directory, qsidd, **overrides).run()

    error = _read_result(directory)["error"]

    assert raised.value.code == code
    assert error["code"] == code

    if code == "qsidd_failed":
        assert (error["exit_code"], error["stderr"]) == (3, "broken")


def test_run_writes_a_failed_result_when_qsidd_is_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("SIST_TRANS_THREE_BIN", str(tmp_path / "missing-qsidd"))

    with pytest.raises(ResultError):
        _runner(tmp_path, FakeQsidd("melting.txt")).run()

    assert _read_result(tmp_path)["error"]["code"] == "qsidd_not_found"


def test_run_passes_the_chosen_text_flags(directory: Path) -> None:
    qsidd = FakeQsidd("melting.txt")

    _runner(
        directory,
        qsidd,
        json_file=None,
        print_parameters=True,
        print_ensemble_average=True,
        nearest_neighbor=True,
        circular=True,
    ).run()

    assert qsidd.commands[0][1:5] == ["-p", "-r", "-n", "-c"]


def test_ir_finder_scorer_returns_the_string_and_version(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class FakeIRFinder:
        """Stand-in for IRFinder that records the shape it was given."""

        irf_version = "3.09"

        def __init__(self, *, temperature: float, shape: str) -> None:
            """Record the settings.

            Args:
                temperature: Temperature in kelvin.
                shape: ``"linear"`` or ``"circular"``.
            """
            self.shape = shape

        def compute_cruciform_energy_string(self, file: str) -> str:
            """Return a string naming the shape and file.

            Args:
                file: Path to the sequence file.

            Returns:
                The string.
            """
            return f"{self.shape},{file}"

    monkeypatch.setattr("sist.runner.IRFinder", FakeIRFinder)
    runner = _runner_without_stand_ins(circular=True)

    repeats = runner._ir_scorer.score("fragment.fa")

    assert repeats == ScoredRepeats("circular,fragment.fa", "3.09")


def test_subprocess_qsidd_runs_a_command() -> None:
    completed = SubprocessQsidd().run([sys.executable, "-c", "print('ran')"])

    assert completed.returncode == 0
    assert completed.stdout == "ran\n"


def _runner_without_stand_ins(**overrides: Any) -> SistRunner:
    """Make a runner with its default IR scorer and qsidd executor.

    Args:
        **overrides: Settings to change from the defaults.

    Returns:
        The runner.
    """
    settings: dict[str, Any] = {
        "file": "fragment.fa",
        "algorithm": "C",
        "temperature": 310.0,
        "superhelical_density": 0.06,
        "salt": 0.01,
        "energy_threshold": 12.0,
        "circular": False,
        "nearest_neighbor": False,
        "print_base_pair": False,
        "print_parameters": False,
        "print_ensemble_average": False,
        "output_file": None,
    }
    settings.update(overrides)

    return SistRunner(**settings)
