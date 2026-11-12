"""Split qsidd's text output into its parts."""

from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass(frozen=True)
class QsiddOutput:
    """The text output of one qsidd run, split into its parts.

    Attributes:
        values: Every ``key = value`` line, with the value kept as text.
        warnings: The text of each ``WARNING:`` line, in order.
        columns: Headings of the per-position table.
        rows: The per-position table, one list of values per position.
    """

    values: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    columns: list[str] = field(default_factory=list)
    rows: list[list[str]] = field(default_factory=list)

    @classmethod
    def parse(cls, text: str) -> QsiddOutput:
        """Split qsidd's text output into values, warnings and the table.

        Most values are printed as ``key = value``, sometimes two to a line
        separated by ``;``. The scaling factor is the one value printed without
        an ``=``.

        Args:
            text: Standard output of a qsidd run.

        Returns:
            The parsed output.
        """
        output = cls()

        for line in text.splitlines():
            if output.columns:
                if line.strip():
                    output.rows.append(line.rstrip("\t").split("\t"))
                continue

            if line.startswith("WARNING: "):
                output.warnings.append(line.removeprefix("WARNING: "))
            elif line.startswith("Position"):
                output.columns.extend(line.split("\t"))
            elif line.startswith("Scaling factor"):
                output.values["Scaling factor"] = line.split()[-1]
            else:
                for part in line.split(";"):
                    if "=" in part:
                        key, value = part.split("=", 1)
                        output.values[key.strip()] = value.strip()

        return output

    def number(self, key: str) -> float | None:
        """Read a printed value as a number.

        Values such as ``-25.4 turns`` keep only the number. ``nan`` and
        ``inf`` are returned as None.

        Args:
            key: Name of the value, as printed before the ``=``.

        Returns:
            The number, or None if it is not a finite number.

        Raises:
            KeyError: If qsidd did not print the value.
        """
        value = float(self.values[key].split()[0])

        return value if math.isfinite(value) else None

    def count(self, key: str) -> int:
        """Read a printed state count as an exact integer.

        Args:
            key: Name of the count, as printed before the ``=``.

        Returns:
            The count.

        Raises:
            KeyError: If qsidd did not print the count.
        """
        return int(float(self.values[key]))

    def column(self, name: str) -> list[str]:
        """Read one column of the per-position table.

        Args:
            name: Column heading, such as ``"P(x)"``.

        Returns:
            The column's values as text, one per position.
        """
        index = self.columns.index(name)

        return [row[index] for row in self.rows]
