"""Unit tests for sist.cli and running sist as a module."""

from __future__ import annotations

import runpy
from typing import Any

import pytest

from sist import cli


class FakeRunner:
    """Stand-in for SistRunner that records how it was made.

    Attributes:
        settings: Keyword arguments of the most recent runner.
        error: Exception to raise when run, if any.
    """

    settings: dict[str, Any] = {}
    error: Exception | None = None

    def __init__(self, **settings: Any) -> None:
        """Record the settings.

        Args:
            **settings: The runner's keyword arguments.
        """
        FakeRunner.settings = settings

    def run(self) -> None:
        """Do nothing, or raise the configured error.

        Raises:
            Exception: The configured error, if set.
        """
        if FakeRunner.error is not None:
            raise FakeRunner.error


@pytest.fixture
def fake_runner(monkeypatch: pytest.MonkeyPatch) -> type[FakeRunner]:
    """Replace SistRunner in the CLI with a stand-in."""
    FakeRunner.settings = {}
    FakeRunner.error = None
    monkeypatch.setattr(cli, "SistRunner", FakeRunner)

    return FakeRunner


def test_main_passes_the_settings_to_the_runner(fake_runner: type[FakeRunner]) -> None:
    cli.main(["-f", "seq.fa", "-a", "Z", "--json", "result.json"])

    assert fake_runner.settings["algorithm"] == "Z"
    assert fake_runner.settings["json_file"] == "result.json"
    assert fake_runner.settings["command"] == [
        "sist",
        "-f",
        "seq.fa",
        "-a",
        "Z",
        "--json",
        "result.json",
    ]


def test_main_exits_on_an_invalid_setting(fake_runner: type[FakeRunner]) -> None:
    with pytest.raises(SystemExit) as raised:
        cli.main(["-f", "seq.fa", "-a", "M", "-T", "0"])

    assert raised.value.code == 1


def test_main_exits_when_the_run_fails(fake_runner: type[FakeRunner]) -> None:
    fake_runner.error = RuntimeError("broken")

    with pytest.raises(SystemExit) as raised:
        cli.main(["-f", "seq.fa", "-a", "M"])

    assert raised.value.code == 1


def test_main_warns_when_run_by_its_old_name(
    fake_runner: type[FakeRunner],
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(cli.sys, "argv", ["master.pl", "-f", "seq.fa", "-a", "M"])

    cli.main()

    assert "deprecated alias" in capsys.readouterr().err
    assert fake_runner.settings["command"][0] == "master.pl"


def test_running_the_package_calls_main(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = []
    monkeypatch.setattr(cli, "main", lambda: calls.append(True))

    runpy.run_module("sist", run_name="__main__")

    assert calls == [True]
