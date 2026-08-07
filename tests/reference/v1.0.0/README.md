# SIST 1.0.0 competition reference

This directory contains reference output for the SIST competition calculation.

## Baseline

* Commit: `dab84fd57f81c3bbb1bb26bd7d4a2fa9674d8a2d`
* Input: `tests/data/pbr322.toy.fa`
* Algorithm: competition (`A`)

## Command

```bash
perl master.pl \
    -f pbr322.toy.fa \
    -a A \
    -o pbr322.toy.compete.txt \
    -b \
    -p \
    -r
```

## Files

* `competition.rebuilt.txt` - output produced from executables rebuilt from the maintained source. This is the regression-test reference.
* `competition.inherited.txt` - output produced by the executable originally included in the repository. This is retained for historical comparison.

All values common to both executables matched exactly, excluding the variable `Run time` line.

The inherited executable also prints 24 additional derived probability fields that are not printed by the rebuilt executable. These are being investigated separately and are not part of the initial regression comparison.

## Regression comparison

The test compares the deterministic calculation results and ignores the `Run time` line.
