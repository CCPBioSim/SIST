"""Record how, where and when a result was made.

A :class:`ProvenanceRecorder` is created when a run starts and finished when its
result is built::

    recorder = ProvenanceRecorder(engine="trans_three", command=["sist", "-a", "M"])
    provenance = recorder.finish(irf_version=None)

``provenance.wall_time_s`` is then the number of seconds between the two.
"""

from __future__ import annotations

import os
import platform
import re
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path

import sist
from sist.results.model import Provenance


class ProvenanceRecorder:
    """Record a run's provenance, starting the clock when it is created.

    The SIST version, git commit, Python version and platform are read when
    :meth:`finish` is called.

    Attributes:
        GIT_SHA_ENV_VAR: Environment variable that sets the git commit, for
            builds made outside a git checkout, such as a conda package.
        engine: qsidd binary the run uses.
        command: The command line as given, or None when run from Python.
        started_at: When the run started.
    """

    GIT_SHA_ENV_VAR = "SIST_GIT_SHA"

    def __init__(self, *, engine: str, command: list[str] | None) -> None:
        """Start the clock for a run, and record its start time.

        Args:
            engine: qsidd binary the run uses, ``"trans_three"`` or
                ``"trans_compete"``.
            command: The command line as given, such as
                ``["sist", "-a", "M", "-f", "seq.fa"]``, or None when run from
                Python.
        """
        self.engine = engine
        self.command = command
        self.started_at = self.timestamp()

        self._start = time.perf_counter()

    def finish(self, *, irf_version: str | None) -> Provenance:
        """Stop the clock and record the provenance.

        Args:
            irf_version: The IRF version, such as ``"3.09"``, for runs that
                score inverted repeats, otherwise None.

        Returns:
            The provenance, with the finish time and the elapsed seconds,
            rounded to milliseconds.
        """
        return Provenance(
            version=sist.__version__,
            git_sha=self.git_sha(),
            engine=self.engine,
            irf_version=irf_version,
            python=platform.python_version(),
            platform=platform.platform(),
            command=self.command,
            started_at=self.started_at,
            finished_at=self.timestamp(),
            wall_time_s=round(time.perf_counter() - self._start, 3),
        )

    @staticmethod
    def timestamp() -> str:
        """The current UTC time in RFC 3339 format.

        Returns:
            The time to the second, such as ``"2026-11-06T10:00:02Z"``.
        """
        return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")

    @classmethod
    def git_sha(cls) -> str | None:
        """Find the git commit SIST is running from.

        Uses the environment variable if it is set. Otherwise reads the commit
        with git when SIST is running from a source checkout, and returns None
        for an installed copy.

        Returns:
            The commit hash, or None if it cannot be found.
        """
        override = os.environ.get(cls.GIT_SHA_ENV_VAR)

        if override:
            return override

        repository_root = Path(__file__).resolve().parents[3]

        if not (repository_root / ".git").exists():
            return None

        try:
            completed = subprocess.run(
                ["git", "-C", str(repository_root), "rev-parse", "HEAD"],
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError:
            return None

        sha = completed.stdout.strip()

        return sha if re.fullmatch(r"[0-9a-f]{7,40}", sha) else None
