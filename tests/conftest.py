from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

SIST_TRANSITIONS = (
    pytest.param(("melting", "M"), id="melting"),
    pytest.param(("z-dna", "Z"), id="z-dna"),
    pytest.param(("cruciform", "C"), id="cruciform"),
)


@dataclass(frozen=True)
class SistRun:
    """Result of running a SIST calculation."""

    name: str
    algorithm: str
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


@pytest.fixture(scope="session")
def built_sist(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """
    Copy and build SIST in a temporary directory.

    This ensures the tests use executables built from the maintained source
    without modifying the repository working tree.
    """

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

    for directory in ("trans_three", "trans_compete"):
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


def run_sist_calculation(
    built_sist: Path,
    *,
    name: str,
    algorithm: str,
) -> SistRun:
    """Run one SIST calculation using the regression test sequence."""

    source_input = built_sist / "tests" / "data" / "pbr322.toy.fa"
    runtime_input = built_sist / "pbr322.toy.fa"

    shutil.copy2(source_input, runtime_input)

    output_directory = built_sist / "test-results"
    output_directory.mkdir(exist_ok=True)

    output_path = output_directory / f"{name}.txt"

    result = run_command(
        [
            "perl",
            "master.pl",
            "-f",
            runtime_input.name,
            "-a",
            algorithm,
            "-o",
            str(output_path.relative_to(built_sist)),
            "-b",
            "-p",
            "-r",
        ],
        cwd=built_sist,
    )

    return SistRun(
        name=name,
        algorithm=algorithm,
        process=result,
        output_path=output_path,
    )


@pytest.fixture(scope="session")
def competition_run(built_sist: Path) -> SistRun:
    """Run the SIST competition calculation once."""

    return run_sist_calculation(
        built_sist,
        name="competition",
        algorithm="A",
    )


@pytest.fixture(
    scope="session",
    params=SIST_TRANSITIONS,
)
def transition_run(
    request: pytest.FixtureRequest,
    built_sist: Path,
) -> SistRun:
    """Run each supported individual SIST transition calculation once."""

    name, algorithm = request.param

    return run_sist_calculation(
        built_sist,
        name=name,
        algorithm=algorithm,
    )
