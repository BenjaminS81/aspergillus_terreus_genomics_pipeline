#!/usr/bin/env bash
set -euo pipefail

QC_ROOT="/mnt/d/Illumina/David_Eisele/terreus_genotyping/terreus_qc_check"
LOG="${QC_ROOT}/complete_data/logs/complete_busco_master.log"
PID_FILE="${QC_ROOT}/complete_data/run_complete_busco.pid"

if [[ -f "${PID_FILE}" ]]; then
    old_pid="$(cat "${PID_FILE}")"
    if kill -0 "${old_pid}" 2>/dev/null; then
        echo "already running with PID ${old_pid}"
        exit 0
    fi
fi

nohup bash "${QC_ROOT}/scripts/run_complete_busco.sh" >> "${LOG}" 2>&1 &
pid="$!"
printf '%s\n' "${pid}" > "${PID_FILE}"
echo "started PID ${pid}"
