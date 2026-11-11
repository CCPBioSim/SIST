"""Unit tests for sist.results.provenance."""

from __future__ import annotations

from datetime import datetime

import pytest

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
