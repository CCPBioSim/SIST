#!/usr/bin/env bash
set -euo pipefail

command -v sist
command -v irf
command -v perl

TEST_ROOT="$(pwd)"

run_reference_test() {
    mode="$1"
    name="$2"
    reference="$3"

    mkdir -p "${name}-test"

    cp \
        "${TEST_ROOT}/tests/data/pbr322.toy.fa" \
        "${name}-test/pbr322.toy.fa"

    (
        cd "${name}-test"

        sist \
            -f pbr322.toy.fa \
            -a "${mode}" \
            -o "${name}.txt" \
            -b \
            -p \
            -r

        grep -v '^Run time =' \
            "${name}.txt" \
            > actual.txt

        grep -v '^Run time =' \
            "${TEST_ROOT}/${reference}" \
            > expected.txt

        diff -u expected.txt actual.txt
    )
}

run_reference_test \
    M \
    melting \
    tests/reference/v1.0.0/melting.txt

run_reference_test \
    Z \
    z-dna \
    tests/reference/v1.0.0/z-dna.txt

run_reference_test \
    C \
    cruciform \
    tests/reference/v1.0.0/cruciform.txt

run_reference_test \
    A \
    competition \
    tests/reference/v1.0.0/competition.rebuilt.txt