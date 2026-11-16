"""Split qsidd's text output into its parts.

:meth:`QsiddOutput.parse` reads the text qsidd prints with ``-b -p -r`` (and
``-P`` for full precision) into named values, warnings and the per-position
table. It does not interpret the values; :mod:`sist.results.reader` does that.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass(frozen=True)
class QsiddOutput:
    """The text output of one qsidd run, split into its parts.

    >>> output = QsiddOutput.parse("Temperature = 310 K\\nPosition\\tP(x)\\n1\\t0.5")
    >>> output.number("Temperature")
    310.0
    >>> output.column("P(x)")
    ['0.5']

    Attributes:
        values: Every ``key = value`` line, keyed by the text before the ``=``
            and with the value kept as text, such as ``{"Temperature": "310 K"}``.
        warnings: The text of each ``WARNING:`` line, after ``WARNING:``, in
            order.
        columns: Headings of the per-position table, such as
            ``["Position", "Base", "P(x)", "G(x)"]``.
        rows: The per-position table, one list of values per position, in the
            same order as ``columns``.
    """

    values: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    columns: list[str] = field(default_factory=list)
    rows: list[list[str]] = field(default_factory=list)

    @classmethod
    def parse(cls, text: str) -> QsiddOutput:
        """Split qsidd's text output into values, warnings and the table.

        Each line is read as one of:

        * ``WARNING: threshold is too small, ...``, added to ``warnings``.
        * ``BAT = 0.258 kcal/mol; BGC = 1.305 kcal/mol``, added to ``values``
          as two entries. Most lines hold one ``key = value``.
        * ``Scaling factor 300``, added to ``values`` without an ``=``.
        * ``Position  Base  P(x)  G(x)``, which starts the per-position table.
          Every non-blank line after it is one row.

        Other lines, such as ``Molecule type: linear DNA``, are skipped.

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

        With ``-P`` qsidd prints counts in full, such as ``1301378``. Without
        it, large counts are printed rounded, such as ``1.30138e+06``, and this
        returns the rounded number.

        Args:
            key: Name of the count, as printed before the ``=``, such as
                ``"Number of three-run states"``.

        Returns:
            The count.

        Raises:
            KeyError: If qsidd did not print the count.
        """
        return int(float(self.values[key]))

    def column(self, name: str) -> list[str]:
        """Read one column of the per-position table.

        Args:
            name: Column heading, such as ``"P(x)"`` or ``"Base"``.

        Returns:
            The column's values as text, one per position, starting at
            position 1.

        Raises:
            ValueError: If the table has no column with that heading.
        """
        index = self.columns.index(name)

        return [row[index] for row in self.rows]
