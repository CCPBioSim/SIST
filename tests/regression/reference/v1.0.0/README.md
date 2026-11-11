# SIST 1.0.0 scientific regression references

This directory contains the reference outputs used to protect the established scientific behaviour of SIST for the 1.0.0 release.

The regression suite covers the individual melting, Z-DNA, and cruciform calculations as well as the competition calculation.

## Baseline

All maintained reference calculations use:

- Input: `tests/data/pbr322.toy.fa`
- Melting algorithm: M
- Z-DNA algorithm: Z
- Cruciform algorithm: C
- Competition algorithm: A

The reference outputs were generated from executables rebuilt from the maintained SIST source.

Repeated calculations produced identical deterministic output. The Run time value is excluded from regression comparison because it varies between executions.

## Commands

These were the commands used to generate the reference outputs below, run
against the `master.pl` pipeline that shipped in SIST 1.0.0 and has since
been replaced by the `sist` Python package. They are kept here as a record
of provenance; the baselines themselves remain the source of truth
regardless of which implementation produced them.

### Melting

    perl master.pl \
        -f pbr322.toy.fa \
        -a M \
        -o melting.txt \
        -b \
        -p \
        -r

### Z-DNA

    perl master.pl \
        -f pbr322.toy.fa \
        -a Z \
        -o z-dna.txt \
        -b \
        -p \
        -r

### Cruciform

    perl master.pl \
        -f pbr322.toy.fa \
        -a C \
        -o cruciform.txt \
        -b \
        -p \
        -r

### Competition

    perl master.pl \
        -f pbr322.toy.fa \
        -a A \
        -o competition.txt \
        -b \
        -p \
        -r

## Files

- melting.txt - maintained reference output for the melting calculation.
- z-dna.txt - maintained reference output for the Z-DNA calculation.
- cruciform.txt - maintained reference output for the cruciform calculation.
- competition.rebuilt.txt - maintained reference output for the competition calculation.
- competition.inherited.txt - output produced by the competition executable originally included in the repository. This is retained for historical provenance and is not the maintained regression reference.

## Scientific regression comparison

Each calculation is run twice on `pbr322.toy.fa`: once with text output and once with a structured JSON result.

The text output must match its reference line for line, with zero tolerance. That covers every reported metric, every sequence position and base, P(x) or P_melt, P_Z and P_cruciform at every position, G(x) where the calculation produces it, and the layout of the output.

The variable Run time line is deliberately ignored because execution time is not scientific output.

The JSON result must hold the same bases, per-position probabilities, G(x) and total state count. The references print 6 significant figures while the JSON holds full precision, so its values are compared to that precision.

## Historical competition output

The maintained competition output and the inherited executable output matched exactly for all values common to both, excluding the variable Run time line.

The inherited executable additionally prints 24 derived probability fields that are not printed by the rebuilt maintained executable. These historical fields are retained for provenance but are not part of the SIST 1.0.0 regression contract.

## Updating references

Reference outputs must not be updated solely to make a failing regression test pass.

Any numerical difference from these baselines should be investigated and understood before a reference is changed.

Intentional scientific changes should document why the established SIST 1.0.0 behaviour has changed.
