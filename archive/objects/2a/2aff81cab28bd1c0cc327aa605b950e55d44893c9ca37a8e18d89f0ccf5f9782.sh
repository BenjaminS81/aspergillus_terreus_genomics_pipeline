#!/usr/bin/env bash
set -euo pipefail
# No job starts without an explicitly selected resource window; do not use while
# the protected Panaroo workload occupies the shared resource budget.
if [[ $# != 2 || "$1" != "--resource-window-open" ]]; then
  echo 'Usage after resource window is established: run_one_bounded.sh --resource-window-open JOB_ID' >&2
  exit 2
fi
TASK_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
TASK_UNIT="terreus-gates-$(date +%s)-$$"
exec systemd-run --user --scope --unit="$TASK_UNIT" \
  -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% \
  /home/benjamin/.local/share/Achromobacter_ANI/envs/busco/bin/python3.12 \
  "$TASK_DIR/run_one_busco.py" "$2" --execute
