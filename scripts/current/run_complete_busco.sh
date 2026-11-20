#!/usr/bin/env bash
set -euo pipefail

QC_ROOT="/mnt/d/Illumina/David_Eisele/terreus_genotyping/terreus_qc_check"
COMPLETE="${QC_ROOT}/complete_data"
BUSCO_ENV="/home/benjamin/.local/share/Achromobacter_ANI/envs/busco"
LINEAGE="${QC_ROOT}/busco_downloads/lineages/eurotiales_odb12.2"
PYTHON="${BUSCO_ENV}/bin/python3.12"
BUSCO="${BUSCO_ENV}/bin/busco"
THREADS_PER_JOB=5
MAX_JOBS=4

mkdir -p "${COMPLETE}/busco_genomes" "${COMPLETE}/busco_transcripts" "${COMPLETE}/logs/genome_busco" "${COMPLETE}/logs/transcript_busco"

run_genome() {
    local label="$1" source="$2" out_dir="${COMPLETE}/busco_genomes/${1}"
    if compgen -G "${out_dir}/short_summary.specific.*.json" > /dev/null; then
        return 0
    fi
    echo "GENOME BUSCO START ${label}"
    local restart=()
    [[ -d "${out_dir}" ]] && restart=(--restart)
    env PATH="${BUSCO_ENV}/bin:/usr/local/bin:/usr/bin:/bin" \
        "${PYTHON}" "${BUSCO}" --in "${source}" --out "${label}" --out_path "${COMPLETE}/busco_genomes" \
        --mode genome --lineage_dataset "${LINEAGE}" --offline --cpu "${THREADS_PER_JOB}" --miniprot \
        --scaffold_composition --tar --opt-out-run-stats "${restart[@]}" \
        > "${COMPLETE}/logs/genome_busco/${label}.log" 2>&1
    echo "GENOME BUSCO DONE ${label}"
}

run_transcript() {
    local label="$1" source="$2" out_dir="${COMPLETE}/busco_transcripts/${1}"
    if compgen -G "${out_dir}/short_summary.specific.*.json" > /dev/null; then
        return 0
    fi
    echo "TRANSCRIPT BUSCO START ${label}"
    local restart=()
    [[ -d "${out_dir}" ]] && restart=(--restart)
    env PATH="${BUSCO_ENV}/bin:/usr/local/bin:/usr/bin:/bin" \
        "${PYTHON}" "${BUSCO}" --in "${source}" --out "${label}" --out_path "${COMPLETE}/busco_transcripts" \
        --mode transcriptome --lineage_dataset "${LINEAGE}" --offline --cpu "${THREADS_PER_JOB}" \
        --tar --opt-out-run-stats "${restart[@]}" \
        > "${COMPLETE}/logs/transcript_busco/${label}.log" 2>&1
    echo "TRANSCRIPT BUSCO DONE ${label}"
}

wait_for_slot() {
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
        wait -n
    done
}

while IFS=$'\t' read -r label species sample data_class has_reads source; do
    [[ "${label}" == "label" ]] && continue
    source="${source%$'\r'}"
    old_json="${QC_ROOT}/busco_results/terreus_genome_busco/${label}.fasta"
    if compgen -G "${old_json}/short_summary.specific.*.json" > /dev/null; then
        continue
    fi
    wait_for_slot
    run_genome "${label}" "${source}" &
done < "${COMPLETE}/complete_genome_manifest.tsv"
wait

while IFS=$'\t' read -r label species sample data_class source matched_genome; do
    [[ "${label}" == "label" ]] && continue
    source="${source%$'\r'}"
    wait_for_slot
    run_transcript "${label}" "${source}" &
done < "${COMPLETE}/complete_transcript_manifest.tsv"
wait

python3 "${QC_ROOT}/scripts/summarize_complete_qc.py" "${QC_ROOT}"
date --iso-8601=seconds > "${COMPLETE}/COMPLETE_DATA_QC_FINISHED.txt"
