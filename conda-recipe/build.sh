#!/usr/bin/env bash
set -euxo pipefail

make -C src/trans_three clean
make -C src/trans_three

make -C src/trans_compete clean
make -C src/trans_compete

install -d "${PREFIX}/libexec/sist/src/trans_three"
install -d "${PREFIX}/libexec/sist/src/trans_compete"

install -m 755 \
    src/trans_three/qsidd \
    "${PREFIX}/libexec/sist/src/trans_three/qsidd"

install -m 755 \
    src/trans_compete/qsidd \
    "${PREFIX}/libexec/sist/src/trans_compete/qsidd"

"${PYTHON}" -m pip install . --no-deps --no-build-isolation -vv
