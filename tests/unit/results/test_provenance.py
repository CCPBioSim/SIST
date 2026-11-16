"""Unit tests for sist.results.provenance."""

from __future__ import annotations

import subprocess
from datetime import datetime
from pathlib import Path

import pytest

from sist.results import provenance
from sist.results.provenance import ProvenanceRecorder


def test_finish_records_the_run(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("SIST_GIT_SHA", "0123456")
    recorder = ProvenanceRecorder(engine="trans_three", command=["sist", "-a", "M"])

    provenance = recorder.finish(irf_version=None)

    assert provenance.engine == "trans_three"
    assert provenance.command == ["sist", "-a", "M"]
    assert provenance.git_sha == "0123456"
    assert provenance.wall_time_s >= 0


def test_timestamp_is_an_rfc_3339_date_time() -> None:
    timestamp = ProvenanceRecorder.timestamp()

    assert timestamp.endswith("Z")
    assert datetime.fromisoformat(timestamp)


def test_to_dict_has_the_schema_shape() -> None:
    provenance = ProvenanceRecorder(engine="trans_compete", command=None).finish(
        irf_version="3.09"
    )

    section = provenance.to_dict()

    assert section["software"]["name"] == "SIST"
    assert section["software"]["irf_version"] == "3.09"
    assert set(section["environment"]) == {"python", "platform"}


def test_git_sha_is_read_from_a_source_checkout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("SIST_GIT_SHA", raising=False)
    module_path = tmp_path / "src" / "sist" / "results" / "provenance.py"
    (tmp_path / ".git").mkdir()
    monkeypatch.setattr(provenance, "__file__", str(module_path))
    monkeypatch.setattr(
        provenance.subprocess,
        "run",
        lambda command, **_: subprocess.CompletedProcess(command, 0, "0123abc\n"),
    )

    assert ProvenanceRecorder.git_sha() == "0123abc"


def test_git_sha_is_none_for_an_installed_copy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("SIST_GIT_SHA", raising=False)
    module_path = tmp_path / "site-packages" / "sist" / "results" / "provenance.py"
    monkeypatch.setattr(provenance, "__file__", str(module_path))

    assert ProvenanceRecorder.git_sha() is None


def test_git_sha_is_none_when_git_cannot_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def missing_git(command: list[str], **_: object) -> None:
        """Behave as if git is not installed.

        Args:
            command: The command that would have run.

        Raises:
            OSError: Always.
        """
        raise OSError("git")

    monkeypatch.delenv("SIST_GIT_SHA", raising=False)
    (tmp_path / ".git").mkdir()
    module_path = tmp_path / "src" / "sist" / "results" / "provenance.py"
    monkeypatch.setattr(provenance, "__file__", str(module_path))
    monkeypatch.setattr(provenance.subprocess, "run", missing_git)

    assert ProvenanceRecorder.git_sha() is None
