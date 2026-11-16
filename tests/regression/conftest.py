from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]

REGRESSION_DATA = Path(__file__).resolve().parent / "data" / "pbr322.toy.fa"


@dataclass(frozen=True)
class RegressionCase:
    """One algorithm run on the regression sequence.

    Attributes:
        name: Name used in test ids and output file names.
        algorithm: The ``-a`` flag.
        reference: File name of its 1.0.0 text reference.
    """

    name: str
    algorithm: str
    reference: str


REGRESSION_CASES = (
    RegressionCase("melting", "M", "melting.txt"),
    RegressionCase("z-dna", "Z", "z-dna.txt"),
    RegressionCase("cruciform", "C", "cruciform.txt"),
    RegressionCase("competition", "A", "competition.rebuilt.txt"),
)


@dataclass(frozen=True)
class SistRun:
    """Result of running a SIST calculation.

    Attributes:
        case: The algorithm that was run.
        process: The finished ``sist`` process.
        output_path: The file it wrote.
    """

    case: RegressionCase
    process: subprocess.CompletedProcess[str]
    output_path: Path


def run_command(
    command: list[str],
    *,
    cwd: Path,
) -> subprocess.CompletedProcess[str]:
    """Run a command and return its completed process."""

    environment = os.environ.copy()
    environment["LC_ALL"] = "C"

    return subprocess.run(
        command,
        cwd=cwd,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def require_success(
    result: subprocess.CompletedProcess[str],
    *,
    description: str,
) -> None:
    """Fail the test setup with useful command output."""

    if result.returncode == 0:
        return

    pytest.fail(
        f"{description} failed with exit code {result.returncode}\n\n"
        f"stdout:\n{result.stdout}\n\n"
        f"stderr:\n{result.stderr}",
        pytrace=False,
    )


def build_sist(
    tmp_path_factory: pytest.TempPathFactory,
) -> Path:
    """Copy and build SIST in a temporary directory."""

    build_root = tmp_path_factory.mktemp("sist-build")
    working_copy = build_root / "SIST"

    shutil.copytree(
        REPOSITORY_ROOT,
        working_copy,
        ignore=shutil.ignore_patterns(
            ".git",
            ".venv",
            ".pytest_cache",
            ".testdata",
            "__pycache__",
            "cmake-build-debug",
        ),
    )

    for directory in ("src/trans_three", "src/trans_compete"):
        clean_result = run_command(
            ["make", "-C", directory, "clean"],
            cwd=working_copy,
        )
        require_success(
            clean_result,
            description=f"Cleaning {directory}",
        )

        build_result = run_command(
            ["make", "-C", directory],
            cwd=working_copy,
        )
        require_success(
            build_result,
            description=f"Building {directory}",
        )

        executable = working_copy / directory / "qsidd"

        if not executable.is_file():
            pytest.fail(
                f"Build completed but {executable} was not created",
                pytrace=False,
            )

    return working_copy


@pytest.fixture(scope="session")
def built_sist_copy(
    tmp_path_factory: pytest.TempPathFactory,
) -> Path:
    """Build SIST's qsidd binaries once, shared by every command fixture."""

    return build_sist(tmp_path_factory)


@pytest.fixture(scope="session")
def sist_command(
    request: pytest.FixtureRequest,
) -> list[str]:
    """
    Return the SIST command under test.

    During conda-build testing, use the installed console script. Otherwise
    run `python -m sist` against the maintained source tree, with its qsidd
    binaries resolved via env vars pointing at a freshly built copy.
    """

    if os.environ.get("CONDA_BUILD_STATE") == "TEST":
        executable = shutil.which("sist")

        if executable is None:
            pytest.fail(
                "The installed sist command was not found in PATH",
                pytrace=False,
            )

        return [executable]

    built_sist_copy: Path = request.getfixturevalue("built_sist_copy")

    os.environ["SIST_TRANS_THREE_BIN"] = str(
        built_sist_copy / "src" / "trans_three" / "qsidd"
    )
    os.environ["SIST_TRANS_COMPETE_BIN"] = str(
        built_sist_copy / "src" / "trans_compete" / "qsidd"
    )

    return [sys.executable, "-m", "sist"]


def run_sist(
    sist_command: list[str],
    tmp_path_factory: pytest.TempPathFactory,
    case: RegressionCase,
    *,
    output_flag: str,
    output_name: str,
) -> SistRun:
    """Run one SIST calculation on the regression sequence.

    Args:
        sist_command: The SIST command under test.
        tmp_path_factory: Factory for the run's working directory.
        case: The algorithm to run.
        output_flag: ``"-o"`` for text output or ``"--json"`` for a result.
        output_name: File name to write.

    Returns:
        The finished run.
    """
    runtime_directory = tmp_path_factory.mktemp(f"sist-{case.name}")
    shutil.copy2(REGRESSION_DATA, runtime_directory)

    process = run_command(
        [
            *sist_command,
            "-f",
            REGRESSION_DATA.name,
            "-a",
            case.algorithm,
            "-b",
            "-p",
            "-r",
            output_flag,
            output_name,
        ],
        cwd=runtime_directory,
    )
    require_success(process, description=f"The {case.name} calculation")

    return SistRun(case, process, runtime_directory / output_name)


@pytest.fixture(scope="session", params=REGRESSION_CASES, ids=lambda case: case.name)
def regression_case(request: pytest.FixtureRequest) -> RegressionCase:
    """Each algorithm, run once per session."""
    case: RegressionCase = request.param

    return case


@pytest.fixture(scope="session")
def text_run(
    regression_case: RegressionCase,
    sist_command: list[str],
    tmp_path_factory: pytest.TempPathFactory,
) -> SistRun:
    """Run an algorithm with text output, as the 1.0.0 references were made."""
    return run_sist(
        sist_command,
        tmp_path_factory,
        regression_case,
        output_flag="-o",
        output_name=f"{regression_case.name}.txt",
    )


@pytest.fixture(scope="session")
def json_run(
    regression_case: RegressionCase,
    sist_command: list[str],
    tmp_path_factory: pytest.TempPathFactory,
) -> SistRun:
    """Run an algorithm with a structured JSON result."""
    return run_sist(
        sist_command,
        tmp_path_factory,
        regression_case,
        output_flag="--json",
        output_name=f"{regression_case.name}.json",
    )
