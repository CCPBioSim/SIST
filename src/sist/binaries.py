"""Locate the compiled qsidd binaries used by sist's orchestration layer."""

from __future__ import annotations

import os
import sys
from pathlib import Path


class BinaryLocator:
    """Locates the compiled qsidd binaries, for either trans_three or trans_compete.

    Attributes:
        ENGINES: For each engine, the environment variable that overrides its
            path and its path inside an installed SIST.
    """

    TRANS_THREE_ENV_VAR = "SIST_TRANS_THREE_BIN"
    TRANS_COMPETE_ENV_VAR = "SIST_TRANS_COMPETE_BIN"

    ENGINES: dict[str, tuple[str, str]] = {
        "trans_three": (TRANS_THREE_ENV_VAR, "libexec/sist/src/trans_three/qsidd"),
        "trans_compete": (
            TRANS_COMPETE_ENV_VAR,
            "libexec/sist/src/trans_compete/qsidd",
        ),
    }

    def binary(self, engine: str) -> Path:
        """Locate the qsidd binary for an engine.

        Args:
            engine: ``"trans_three"`` or ``"trans_compete"``.

        Returns:
            Path to the qsidd executable.

        Raises:
            KeyError: If the engine is not known.
        """
        env_var, installed_relative_path = self.ENGINES[engine]

        return self._resolve(env_var, installed_relative_path)

    def trans_three_binary(self) -> Path:
        """Locate the melting/Z-DNA/cruciform qsidd binary.

        Returns:
            Path to the trans_three qsidd executable.
        """
        return self.binary("trans_three")

    def trans_compete_binary(self) -> Path:
        """Locate the competition qsidd binary.

        Returns:
            Path to the trans_compete qsidd executable.
        """
        return self.binary("trans_compete")

    @staticmethod
    def _resolve(env_var: str, installed_relative_path: str) -> Path:
        """Resolve a qsidd binary path via an env var override or install layout.

        Args:
            env_var: Environment variable that, if set, is used directly as
                the binary path.
            installed_relative_path: Path to the binary relative to
                ``sys.prefix``, matching the conda package's install layout.

        Returns:
            Path to an existing qsidd executable.

        Raises:
            FileNotFoundError: If neither the env var nor the install
                layout resolves to an existing file.
        """
        override = os.environ.get(env_var)
        path = (
            Path(override) if override else Path(sys.prefix) / installed_relative_path
        )

        if not path.is_file():
            raise FileNotFoundError(
                f"Could not find a qsidd executable at {path}. Set {env_var} to "
                "its path, or install sist so it is available there."
            )

        return path
