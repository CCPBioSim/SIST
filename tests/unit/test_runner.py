"""Unit tests for sist.runner, using stand-ins for IRF and qsidd."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest

from sist.results.diagnostics import ErrorCode, ResultError
from sist.runner import ScoredRepeats, SistRunner
from sist.schema import ResultValidator

QSIDD_PATH = Path(__file__).resolve().parent / "data" / "qsidd"


class FakeQsidd:
    """Stand-in for qsidd that prints a saved output and records each command.

    Attributes:
        output_name: File name of the saved qsidd output to print.
        returncode: Exit code to report.
        stderr: Standard error to report.
        commands: Every command line it was asked to run.
    """

    def __init__(self, output_name: str, *, returncode: int = 0, stderr: str = ""):
        """Set up the stand-in.

        Args:
            output_name: File name of the saved qsidd output to print.
            returncode: Exit code to report.
            stderr: Standard error to report.
        """
        self.output_name = output_name
        self.returncode = returncode
        self.stderr = stderr
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
def binaries(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point SIST at placeholder qsidd binaries, so it can find them."""
    placeholder = tmp_path / "qsidd"
    placeholder.write_text("", encoding="utf-8")
    monkeypatch.setenv("SIST_TRANS_THREE_BIN", str(placeholder))
    monkeypatch.setenv("SIST_TRANS_COMPETE_BIN", str(placeholder))

    return tmp_path


def _runner(directory: Path, qsidd: FakeQsidd, **overrides: Any) -> SistRunner:
    """Make a runner for the 100 bp fragment, using the stand-ins.

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


def test_run_writes_a_valid_result(binaries: Path) -> None:
    _runner(binaries, FakeQsidd("melting.txt")).run()

    result = _read_result(binaries)

    assert ResultValidator().validate(result) == []
    assert result["provenance"]["command"] == ["sist", "-a", "M", "-f", "fragment.fa"]


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
    binaries: Path, flag: str, output_name: str, algorithm_flags: list[str]
) -> None:
    qsidd = FakeQsidd(output_name)

    _runner(binaries, qsidd, algorithm=flag).run()

    command = qsidd.commands[0]

    assert command[1:5] == ["-b", "-p", "-r", "-P"]
    assert [part for part in command if part in ("-Z", "-C", "-X")] == algorithm_flags
    assert ResultValidator().validate(_read_result(binaries)) == []


def test_run_passes_nearest_neighbor_only_where_supported(binaries: Path) -> None:
    qsidd = FakeQsidd("z-dna.txt")

    _runner(binaries, qsidd, algorithm="Z", nearest_neighbor=True).run()

    assert "-n" not in qsidd.commands[0]


def test_run_with_only_json_prints_no_text(
    binaries: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    _runner(binaries, FakeQsidd("melting.txt")).run()

    assert capsys.readouterr().out == ""


def test_run_writes_text_as_well_when_an_output_file_is_given(binaries: Path) -> None:
    qsidd = FakeQsidd("melting.txt")
    output_path = binaries / "result.txt"

    _runner(binaries, qsidd, output_file=str(output_path), print_base_pair=True).run()

    assert output_path.read_text(encoding="utf-8").startswith("DNA sequence file")
    assert "-P" in qsidd.commands[0]
    assert "-P" not in qsidd.commands[1]
    assert "-b" in qsidd.commands[1]


def test_run_without_json_prints_the_text_output(
    binaries: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    _runner(binaries, FakeQsidd("melting.txt"), json_file=None).run()

    expected = (QSIDD_PATH / "melting.txt").read_text(encoding="utf-8")

    assert capsys.readouterr().out == expected
    assert not (binaries / "result.json").exists()


def test_run_without_json_raises_when_qsidd_fails(binaries: Path) -> None:
    runner = _runner(binaries, FakeQsidd("melting.txt", returncode=3), json_file=None)

    with pytest.raises(subprocess.CalledProcessError):
        runner.run()


def test_run_rejects_an_unknown_algorithm(binaries: Path) -> None:
    with pytest.raises(ValueError):
        _runner(binaries, FakeQsidd("melting.txt"), algorithm="X").run()


@pytest.mark.parametrize(
    ("overrides", "qsidd", "code"),
    [
        ({}, FakeQsidd("melting.txt", returncode=3, stderr="broken"), "qsidd_failed"),
        ({}, FakeQsidd("fragment.fa"), "internal_error"),
        (
            {"algorithm": "C", "ir_scorer": FakeScorer(missing_irf=True)},
            FakeQsidd("cruciform.txt"),
            "irf_not_found",
        ),
    ],
)
def test_run_writes_a_failed_result(
    binaries: Path, overrides: dict[str, Any], qsidd: FakeQsidd, code: str
) -> None:
    with pytest.raises(ResultError) as raised:
        _runner(binaries, qsidd, **overrides).run()

    result = _read_result(binaries)

    assert raised.value.code == code
    assert ResultValidator().validate(result) == []
    assert result["error"]["code"] == code


def test_run_records_qsidd_exit_code_and_stderr(binaries: Path) -> None:
    qsidd = FakeQsidd("melting.txt", returncode=3, stderr="broken")

    with pytest.raises(ResultError):
        _runner(binaries, qsidd).run()

    error = _read_result(binaries)["error"]

    assert error["exit_code"] == 3
    assert error["stderr"] == "broken"


def test_run_writes_a_failed_result_when_qsidd_is_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("SIST_TRANS_THREE_BIN", str(tmp_path / "missing-qsidd"))

    with pytest.raises(ResultError):
        _runner(tmp_path, FakeQsidd("melting.txt")).run()

    assert _read_result(tmp_path)["error"]["code"] == ErrorCode.QSIDD_NOT_FOUND


def test_run_writes_a_failed_result_for_a_missing_sequence_file(
    binaries: Path,
) -> None:
    runner = _runner(
        binaries, FakeQsidd("melting.txt"), file=str(binaries / "missing.fa")
    )

    with pytest.raises(ResultError):
        runner.run()

    assert _read_result(binaries)["error"]["code"] == ErrorCode.INPUT_NOT_FOUND


def test_run_writes_an_internal_error_when_qsidd_prints_nan(binaries: Path) -> None:
    """A NaN where a number is required is reported, never written as a result."""
    text = (QSIDD_PATH / "melting.txt").read_text(encoding="utf-8")
    broken = binaries / "broken.txt"
    broken.write_text(
        text.replace("Average number of runs = ", "Average number of runs = -nan #"),
        encoding="utf-8",
    )
    qsidd = FakeQsidd(str(broken.resolve()))

    with pytest.raises(ResultError) as raised:
        _runner(binaries, qsidd).run()

    assert raised.value.code == ErrorCode.INTERNAL_ERROR
    assert "is not of type 'number'" in _read_result(binaries)["error"]["message"]
