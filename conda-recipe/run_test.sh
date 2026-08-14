#!/usr/bin/env bash
set -euo pipefail

command -v sist
command -v irf
command -v perl
command -v python
command -v pytest

python -m pytest tests -vv