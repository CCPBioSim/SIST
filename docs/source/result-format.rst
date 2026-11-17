Structured Results
==================

This page describes version 1 of the SIST structured result: what a SIST
result contains, and how it is written as JSON.

Writing a result
----------------

Add ``--json`` to any ``sist`` command to write the result to a JSON file:

.. code-block:: bash

   sist -a M -f sequence.fa --json result.json

With ``--json``, the text output is only written when ``-o`` is also given:

.. code-block:: bash

   sist -a M -f sequence.fa --json result.json -o result.txt

A JSON result always holds every value SIST works out, whatever ``-b``,
``-p`` and ``-r`` are set to. qsidd is run at full precision for it, so the
numbers can have more digits than the text output.

If the run fails, a result is still written, with ``"status": "failed"`` and
the reason in ``error``. See `Failed runs`_.

Schema
------

The JSON schema is :download:`sist-result-v1.schema.json
<../../src/sist/schemas/sist-result-v1.schema.json>`. It ships inside the
``sist`` package as ``sist/schemas/sist-result-v1.schema.json``.

Validating a result
-------------------

``sist`` checks every result before writing it. To check a result file
yourself, from the command line:

.. code-block:: bash

   sist-validate result.json

``sist-validate`` takes any number of files, prints ``valid`` or the problems
for each one, and exits with 1 if any file is not valid.

From Python:

.. code-block:: python

   import json
   from sist.schema import ResultValidator

   with open("result.json", encoding="utf-8") as file:
       result = json.load(file)

   for issue in ResultValidator().validate(result):
       print(issue)

``ResultValidator.validate`` returns a list of problems, which is empty if the
result is valid. ``ResultValidator.check`` raises ``ResultValidationError``
instead.

Validation happens in two parts:

1. **The JSON schema** checks the shape of the result: types, ranges, required
   keys and allowed names. Any JSON Schema tool can run it.
2. **Cross-checks** in ``ResultValidator`` compare one part of a result with
   another, which JSON Schema cannot express:

   * ``calculation.transitions`` matches the algorithm.
   * The keys of ``transitions`` match ``calculation.transitions``.
   * ``summary.pairwise`` only appears with more than one transition, and
     never pairs a transition with itself.
   * ``input.analysed_sequence`` and every per-position array are as long as
     the sequence.
   * ``summary.state_counts.by_number_of_runs`` adds up to ``total``.
   * ``ln_partition_function`` matches ``partition_function``.
   * The timestamps are real dates.

   Tools that only run the schema miss these, so use ``ResultValidator`` or
   ``sist-validate`` where you can.

Layout
------

A result is one JSON object. These keys are always present:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Contents
   * - ``schema_version``
     - Version of this model, for example ``"1.0.0"``.
   * - ``status``
     - ``"success"`` or ``"failed"``.
   * - ``provenance``
     - How, where and when the result was made. See `Provenance`_.
   * - ``calculation``
     - The ``algorithm`` (M, Z, C or A) and the ``transitions`` it calculates.
   * - ``input``
     - Sequence file, FASTA header, length, the analysed sequence and its
       SHA-256, and molecule type.
   * - ``parameters``
     - Temperature, superhelical density, salt, energy threshold and melting
       energetics, as used by the calculation. C and A runs also have the salt
       the IR finder used.
   * - ``derived_parameters``
     - Linking difference, torsional stiffness, and the copolymeric A-T and
       G-C melting energies.
   * - ``summary``
     - Whole-molecule results: partition function, state counts and average
       number of runs. Calculations with more than one transition also have
       ``pairwise`` probabilities.
   * - ``transitions``
     - The results for each transition. See `Transitions`_.
   * - ``warnings``
     - List of warnings. Empty if there are none.
   * - ``error``
     - Why the run failed, or ``null``.

These keys are optional:

* ``inverted_repeats``: the cruciform candidates scored by the IR finder
  (algorithms C and A), each with its ``start``, ``length_bp`` and
  ``energy_kcal_per_mol``.
* ``diagnostics``: values for checking a run, such as the internal scaling
  factor and the energy window. See `Versioning`_.

Every field has a description in the schema, so editors that understand JSON
Schema show it as hover text.

Conventions
-----------

Transitions
^^^^^^^^^^^

Transitions are named ``melting``, ``z_dna`` and ``cruciform``. Each one
calculated has an entry in ``transitions``, keyed by its name:

.. code-block:: text

   "transitions": {
     "melting": {
       "nucleation_energy_kcal_per_mol": 10.84,
       "probability": 0.9999999999924621,
       "expected_runs": 1.9404831456507936,
       "expected_transformed_bp": 79.56442046580413,
       "profile": {
         "probability": [2.439586212219559e-07, ...],
         "free_energy_kcal_per_mol": [11.056580176644488, ...]
       }
     }
   }

``probability`` is the chance of at least one run of the transition.
``profile.probability`` is P(x), the chance that each position is in the
transition. ``profile.free_energy_kcal_per_mol`` is G(x), and is only present
for melting runs (algorithm M).

Single-transition and competition results have the same layout, so the same
code reads both.

Pairwise probabilities
^^^^^^^^^^^^^^^^^^^^^^

Competition results have ``summary.pairwise``, with three matrices keyed by
transition name. Each is read as ``matrix[a][b]``:

* ``joint_probability``: P(a and b).
* ``conditional_probability``: P(a given b).
* ``conditional_probability_given_absent``: P(a given not b).

A transition is never paired with itself.

Positions and arrays
^^^^^^^^^^^^^^^^^^^^

* Sequence positions are 1-based and inclusive, as in the text output and
  IRF.
* Every per-position array has exactly ``input.sequence_length_bp`` entries.
  Index ``k`` (0-based, as in JSON and Python) holds position ``k + 1``.
* ``input.analysed_sequence`` is a string with one letter per position.
* ``summary.state_counts.by_number_of_runs`` uses the same idea: entry ``k``
  counts states with ``k + 1`` runs.
* An inverted repeat covers positions ``start`` to ``start + length_bp - 1``.
  For circular molecules it can wrap past the end of the sequence.
* Positions refer to the analysed sequence. qsidd drops characters other
  than A, C, G, T and N, which shifts every later position. SIST warns when
  this happens.

Units
^^^^^

Units are part of the key name wherever a value has one.

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Quantity
     - Unit
     - Key suffix
   * - Temperature
     - K
     - ``_k``
   * - Salt concentration
     - mol/L
     - ``_m``
   * - Energies, including G(x) and the threshold
     - kcal/mol
     - ``_kcal_per_mol``
   * - Torsional stiffness
     - kcal bp / rad\ :sup:`2`
     - ``_kcal_bp_per_rad2``
   * - Linking difference
     - turns
     - ``_turns``
   * - Lengths and expected base pair counts
     - base pairs
     - ``_bp``
   * - Wall time
     - seconds
     - ``_s``
   * - Probabilities, superhelical density, partition function, run counts
     - dimensionless
     - none

Numbers
^^^^^^^

* Floating point values are written with enough digits to read back the
  exact double.
* State counts are exact integers. The text output prints large counts in
  scientific notation, such as ``2.87355e+06``, which loses digits.
* ``ln_partition_function`` is always present in a successful result.
  ``partition_function`` is ``null`` when it is too small to store as a
  double, which happens at high superhelical density.

Missing values
^^^^^^^^^^^^^^

* A value that is undefined, or that does not apply, is ``null``. For
  example, P(a given b) is ``null`` when b never happens, and
  ``bat_kcal_per_mol`` is ``null`` for a Z-DNA run.
* Required keys are always written, even when their value is ``null``. This
  separates "does not apply" from "forgot to write it".
* ``NaN`` and ``Infinity`` are never written. They are not valid JSON.
* G(x) is ``null`` where the melting probability is 0. Elsewhere it is the
  value qsidd prints, which qsidd caps at 12 kcal/mol.

Provenance
----------

``provenance`` records how, where and when a result was made, as opposed to
what was calculated:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Field
     - Contents
   * - ``software``
     - ``name``, ``version``, ``git_sha`` (null for installed packages), the
       qsidd ``engine`` used, and ``irf_version`` for C and A runs.
   * - ``environment``
     - Python version and platform.
   * - ``command``
     - The command line as given, one argument per entry. ``null`` when SIST
       is used from Python.
   * - ``started_at``, ``finished_at``, ``wall_time_s``
     - UTC start and finish times, and elapsed seconds.

Everything in ``provenance`` can change between machines or runs without the
science changing. Comparisons between results, including regression tests,
should drop the whole block and compare the rest exactly.

Warnings
--------

Each warning has a ``code``, a ``message``, a ``source`` (``sist``,
``ir_finder`` or ``qsidd``) and optional ``details``. Match on ``code``, as
the message can change. New codes can appear in any minor version, so readers
must accept codes they do not know.

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Code
     - Source
     - Raised when
   * - ``threshold_too_low``
     - qsidd
     - The energy threshold is below 9 kcal/mol.
   * - ``threshold_too_high``
     - qsidd
     - The energy threshold is above 15 kcal/mol.
   * - ``superhelical_density_out_of_range``
     - qsidd
     - The absolute superhelical density is above 0.15.
   * - ``temperature_out_of_range``
     - qsidd
     - The temperature is below 220 K or above 320 K.
   * - ``salt_concentration_out_of_range``
     - qsidd
     - The salt concentration is below 0.0001 M.
   * - ``sequence_too_short``
     - sist
     - The sequence is shorter than 1,500 bp.
   * - ``sequence_too_long``
     - sist
     - The sequence is longer than 10,000 bp.
   * - ``ambiguous_bases_replaced``
     - sist
     - The sequence has ``N`` bases, which are analysed as ``G``.
   * - ``unrecognised_characters_ignored``
     - sist
     - The sequence has characters other than A, C, G, T and N. They are
       dropped, which shifts later positions.
   * - ``nearest_neighbor_ignored``
     - sist
     - ``-n`` was given for a Z-DNA or cruciform run.
   * - ``cruciform_salt_fixed``
     - sist
     - The salt is not 0.01 M in a C or A run. The IR finder still scores
       inverted repeats at 0.01 M.
   * - ``no_inverted_repeats_found``
     - ir_finder
     - IRF found no candidate inverted repeats.

qsidd checks the sequence length against the size of the file, including the
header and newlines, rather than the number of base pairs. So SIST raises the
two length warnings itself, from the base pair count.

Failed runs
-----------

A failed run still writes a result, so a script can tell a failure apart from
a missing file:

* ``status`` is ``"failed"``.
* ``error`` has a ``code``, a ``message`` and the ``stage`` that failed. It can
  also hold the external program's ``exit_code`` and ``stderr``, and
  ``details``.
* ``summary`` and ``transitions`` are ``null``.
* The input, parameters and derived parameters are filled in as far as the
  run got. Anything not known is ``null``.
* Warnings raised before the failure are kept.
* ``sist`` exits with a non-zero code.

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Code
     - Stage
     - Raised when
   * - ``input_not_found``
     - input
     - The sequence file does not exist.
   * - ``empty_sequence``
     - input
     - The sequence file has no bases.
   * - ``multiple_sequences``
     - input
     - The file has more than one FASTA record. qsidd would join them into
       one sequence while the IR finder reads only the first.
   * - ``irf_not_found``
     - ir_finder
     - IRF is not installed, and the algorithm is C or A.
   * - ``qsidd_not_found``
     - qsidd
     - The qsidd binary cannot be found.
   * - ``qsidd_failed``
     - qsidd
     - qsidd exited with an error.
   * - ``partition_function_underflow``
     - qsidd
     - The partition function underflowed to zero. qsidd still exits normally,
       but every probability is meaningless.
   * - ``internal_error``
     - internal
     - Anything else, including a result that fails validation.

The schema also allows the stages ``validation`` and ``output``, for codes
added in later minor versions.

There is no partial status. A run either produces all of its results or is
marked as failed.

Versioning
----------

``schema_version`` follows ``MAJOR.MINOR.PATCH``:

* **Major**: a change that can break a reader. Removing or renaming a field,
  changing its type, unit, meaning or indexing, or making an optional field
  required.
* **Minor**: a change existing readers can ignore. Adding an optional field,
  a transition, an algorithm, a warning code or an error code.
* **Patch**: wording and documentation only.

Each major version has one schema file, named after it
(``sist-result-v1.schema.json``), with a matching ``$id``. The file is updated
in place for minor and patch releases, and accepts any ``1.x.y`` version. A
new major version gets a new file, and the old one stays available.

Readers should check that the major version is one they support, and ignore
fields they do not know. Writers validate against the latest schema for their
major version, which rejects unknown fields to catch mistakes.

``diagnostics`` is outside these rules. Its fields can change in any release,
so do not rely on them.

Text output mapping
-------------------

Where each value in the text output goes:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Text output
     - Result field
   * - ``DNA sequence file``
     - ``input.sequence_file``
   * - ``Sequence Length``
     - ``input.sequence_length_bp``
   * - ``Transition: ...`` or ``... melting energetics``
     - ``calculation.transitions``, ``parameters.melting_energetics``
   * - ``Molecule type``
     - ``input.molecule_type``
   * - ``Stress Level``
     - ``parameters.superhelical_density``
   * - ``linking difference``
     - ``derived_parameters.linking_difference_turns``
   * - ``Threshold``, ``Temperature``, ``Salt Concentration``
     - ``parameters.energy_threshold_kcal_per_mol``, ``temperature_k``,
       ``salt_concentration_m``
   * - ``K`` or ``C``
     - ``derived_parameters.torsional_stiffness_kcal_bp_per_rad2``
   * - ``Nucleation energy`` (one per transition in competition)
     - ``transitions.<name>.nucleation_energy_kcal_per_mol``
   * - ``BAT``, ``BGC``
     - ``derived_parameters.bat_kcal_per_mol``, ``bgc_kcal_per_mol``
   * - ``MinE``, ``Max_E``, ``Scaling factor``
     - ``diagnostics.energy_window``, ``diagnostics.scaling_factor``
   * - ``Number of one-run states`` ... ``four-run states``
     - ``summary.state_counts.by_number_of_runs``
   * - ``Total number of states``
     - ``summary.state_counts.total``
   * - ``Total Partition Function``
     - ``summary.partition_function``, ``summary.ln_partition_function``
   * - ``Average number of runs``
     - ``summary.average_number_of_runs``
   * - ``Transition Probability``, ``Prob_M``, ``Prob_Z``, ``Prob_C``
     - ``transitions.<name>.probability``
   * - ``Number of M_runs`` (and Z, C)
     - ``transitions.<name>.expected_runs``
   * - ``Number of melted bases``, ``N_Melting`` (and Z, C)
     - ``transitions.<name>.expected_transformed_bp``
   * - ``Run time``
     - ``provenance.wall_time_s``
   * - ``Base`` column
     - ``input.analysed_sequence``
   * - ``P(x)``, ``P_melt``, ``P_Z``, ``P_cruciform`` columns
     - ``transitions.<name>.profile.probability``
   * - ``G(x)`` column
     - ``transitions.melting.profile.free_energy_kcal_per_mol``
   * - ``WARNING: ...`` lines
     - ``warnings``

Values SIST works out but does not print in the text output:

* Joint and conditional probabilities for pairs of transitions in competition
  runs, in ``summary.pairwise``.
* The scored inverted repeats passed to qsidd, in ``inverted_repeats``.
* The salt used by the IR finder, in
  ``parameters.cruciform_salt_concentration_m``.
