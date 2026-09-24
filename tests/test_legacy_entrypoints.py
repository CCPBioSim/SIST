"""Tests for the deprecated `master.pl`/`IR_finder.pl` console-script aliases.

These exist for scripted workflows that still invoke the old Perl script
names directly; both aliases will be removed in the next release.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest
from conftest import REPOSITORY_ROOT

pytestmark = pytest.mark.regression


def test_master_pl_alias_matches_sist_output(
    request: pytest.FixtureRequest,
    tmp_path_factory: pytest.TempPathFactory,
) -> None:
    """The deprecated `master.pl` alias should behave exactly like `sist`.

    `built_sist_copy` is looked up lazily (rather than taken as a normal
    fixture parameter) so conda-build testing never triggers it: that build
    only has `tests/` and `pyproject.toml` available, not the C++ sources.
    """

    environment = os.environ.copy()

    if os.environ.get("CONDA_BUILD_STATE") != "TEST":
        built_sist_copy: Path = request.getfixturevalue("built_sist_copy")
        environment["SIST_TRANS_THREE_BIN"] = str(
            built_sist_copy / "src" / "trans_three" / "qsidd"
        )
        environment["SIST_TRANS_COMPETE_BIN"] = str(
            built_sist_copy / "src" / "trans_compete" / "qsidd"
        )

    runtime_directory = tmp_path_factory.mktemp("master-pl-alias")
    input_path = runtime_directory / "pbr322.toy.fa"
    shutil.copy2(REPOSITORY_ROOT / "tests" / "data" / "pbr322.toy.fa", input_path)

    result = subprocess.run(
        ["master.pl", "-f", input_path.name, "-a", "M", "-b", "-p", "-r"],
        cwd=runtime_directory,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "deprecated" in result.stderr.lower()
    assert "Sequence Length" in result.stdout


@pytest.mark.skipif(shutil.which("irf") is None, reason="requires the irf executable")
def test_ir_finder_pl_alias_keeps_original_interface(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The deprecated `IR_finder.pl` alias should keep its positional-arg interface."""

    monkeypatch.chdir(tmp_path)
    Path("toy.fa").write_text(
        ">toy.fa\nAAACCACCGCTTTTTTTTTTTTTTTTTTTTTTTTTGGCGGTGGTTT\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        ["IR_finder.pl", "310.0", "linear", "toy.fa"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "deprecated" in result.stderr.lower()
    assert result.stdout.strip().startswith("1,0,10000,|")
