#!/bin/sh
# Download the public tau-bench historical trajectories (MIT license) at a pinned
# commit. Read-only download from GitHub; nothing is uploaded.
set -eu
SHA=59a200c6d575d595120f1cb70fea53cef0632f6b
BASE="https://raw.githubusercontent.com/sierra-research/tau-bench/$SHA"
mkdir -p data/raw/tau_bench
for f in gpt-4o-airline.json gpt-4o-retail.json sonnet-35-new-airline.json sonnet-35-new-retail.json; do
  [ -f "data/raw/tau_bench/$f" ] || curl -sSfL -o "data/raw/tau_bench/$f" "$BASE/historical_trajectories/$f"
done
[ -f data/raw/tau_bench/LICENSE ] || curl -sSfL -o data/raw/tau_bench/LICENSE "$BASE/LICENSE"
echo "tau-bench files in data/raw/tau_bench (commit $SHA)"
