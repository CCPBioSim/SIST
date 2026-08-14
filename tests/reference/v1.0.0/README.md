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

The regression suite compares deterministic scientific output against these references with zero relative and absolute tolerance.

For individual transition calculations, the comparison includes:

- Reported deterministic scientific metrics.
- Every sequence position.
- The base reported at every position.
- P(x) at every position.
- G(x) where it is produced by the calculation.

For the competition calculation, the comparison includes:

- Reported deterministic scientific metrics.
- Every sequence position.
- The base reported at every position.
- P_melt, P_Z, and P_cruciform at every position.

The structure of the scientific profile output is also checked so that unexpected changes to the reported columns cause the regression tests to fail.

The variable Run time line is deliberately ignored because execution time is not scientific output.

## Historical competition output

The maintained competition output and the inherited executable output matched exactly for all values common to both, excluding the variable Run time line.

The inherited executable additionally prints 24 derived probability fields that are not printed by the rebuilt maintained executable. These historical fields are retained for provenance but are not part of the SIST 1.0.0 regression contract.

## Updating references

Reference outputs must not be updated solely to make a failing regression test pass.

Any numerical difference from these baselines should be investigated and understood before a reference is changed.

Intentional scientific changes should document why the established SIST 1.0.0 behaviour has changed.