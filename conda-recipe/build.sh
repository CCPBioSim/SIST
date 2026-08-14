#!/usr/bin/env bash
set -euxo pipefail

make -C trans_three clean
make -C trans_three

make -C trans_compete clean
make -C trans_compete

install -d "${PREFIX}/bin"
install -d "${PREFIX}/libexec/sist"
install -d "${PREFIX}/libexec/sist/trans_three"
install -d "${PREFIX}/libexec/sist/trans_compete"

install -m 755 \
    master.pl \
    "${PREFIX}/libexec/sist/master.pl"

install -m 755 \
    IR_finder.pl \
    "${PREFIX}/libexec/sist/IR_finder.pl"

install -m 755 \
    trans_three/qsidd \
    "${PREFIX}/libexec/sist/trans_three/qsidd"

install -m 755 \
    trans_compete/qsidd \
    "${PREFIX}/libexec/sist/trans_compete/qsidd"

cat > "${PREFIX}/bin/sist" << 'EOF'
#!/usr/bin/env bash

PREFIX="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

exec "${PREFIX}/bin/perl" \
    "${PREFIX}/libexec/sist/master.pl" \
    "$@"
EOF

chmod 755 "${PREFIX}/bin/sist"