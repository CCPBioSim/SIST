"""Locate the compiled qsidd binaries used by sist's orchestration layer."""

from __future__ import annotations

import os
import sys
from pathlib import Path


class BinaryLocator:
    """Locates the compiled qsidd binaries, for either trans_three or trans_compete."""

    TRANS_THREE_ENV_VAR = "SIST_TRANS_THREE_BIN"
    TRANS_COMPETE_ENV_VAR = "SIST_TRANS_COMPETE_BIN"

    def trans_three_binary(self) -> Path:
        """Locate the melting/Z-DNA/cruciform qsidd binary.

        Returns:
            Path to the trans_three qsidd executable.
        """
        return self._resolve(
            self.TRANS_THREE_ENV_VAR, "libexec/sist/src/trans_three/qsidd"
        )

    def trans_compete_binary(self) -> Path:
        """Locate the competition qsidd binary.

        Returns:
            Path to the trans_compete qsidd executable.
        """
        return self._resolve(
            self.TRANS_COMPETE_ENV_VAR, "libexec/sist/src/trans_compete/qsidd"
        )

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
