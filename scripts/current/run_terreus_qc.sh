#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/d/Illumina/David_Eisele/terreus_genotyping"
ASSEMBLY_ROOT="${PROJECT_ROOT}/results/0_1_spades"
OUT_ROOT="${PROJECT_ROOT}/terreus_qc_check"
BUSCO_ENV="/home/benjamin/.local/share/Achromobacter_ANI/envs/busco"
SEQKIT="/home/benjamin/miniforge3/envs/klebsiella_pangenome/bin/seqkit"
MINIMAP2="/home/benjamin/miniforge3/envs/klebsiella_analysis/bin/minimap2"
SAMTOOLS="/home/benjamin/miniforge3/envs/klebsiella_analysis/bin/samtools"
READ_ROOT="${PROJECT_ROOT}/results/0_0_fastp"

mkdir -p \
    "${OUT_ROOT}/logs" \
    "${OUT_ROOT}/metrics" \
    "${OUT_ROOT}/read_mapping/per_sample" \
    "${OUT_ROOT}/busco_downloads" \
    "${OUT_ROOT}/busco_results"

MANIFEST="${OUT_ROOT}/assembly_manifest.tsv"
printf 'label\tspecies_group\tsample\tsource_fasta\n' > "${MANIFEST}"

STAGING="$(mktemp -d /tmp/terreus_qc_inputs.XXXXXX)"
cleanup() {
    rm -rf -- "${STAGING}"
}
trap cleanup EXIT

while IFS= read -r assembly; do
    sample="$(basename "$(dirname "${assembly}")")"
    species="$(basename "$(dirname "$(dirname "${assembly}")")")"
    label="${species}__${sample}"
    printf '%s\t%s\t%s\t%s\n' "${label}" "${species}" "${sample}" "${assembly}" >> "${MANIFEST}"
    ln -s "${assembly}" "${STAGING}/${label}.fasta"
done < <(find "${ASSEMBLY_ROOT}" -mindepth 3 -maxdepth 3 -type f -name 'scaffolds.fasta' -print | sort)

ASSEMBLY_COUNT="$(( $(wc -l < "${MANIFEST}") - 1 ))"
if [[ "${ASSEMBLY_COUNT}" -ne 28 ]]; then
    echo "Expected 28 final scaffolds.fasta files but found ${ASSEMBLY_COUNT}" >&2
    exit 2
fi

python3 "${OUT_ROOT}/scripts/assembly_metrics.py" "${MANIFEST}" "${OUT_ROOT}/metrics/assembly_metrics.tsv"

"${SEQKIT}" stats --all --tabular "${STAGING}"/*.fasta > "${OUT_ROOT}/metrics/seqkit_stats.tsv" 2> "${OUT_ROOT}/logs/seqkit.log"

READ_MANIFEST="${OUT_ROOT}/read_mapping/read_manifest.tsv"
printf 'label\tassembly\tread1\tread2\n' > "${READ_MANIFEST}"
while IFS=$'\t' read -r label species sample assembly; do
    [[ "${label}" == "label" ]] && continue
    read1="${READ_ROOT}/${species}/${sample}_1.fastq.gz"
    read2="${READ_ROOT}/${species}/${sample}_2.fastq.gz"
    if [[ ! -s "${read1}" || ! -s "${read2}" ]]; then
        echo "Missing trimmed read pair for ${label}" >&2
        exit 3
    fi
    printf '%s\t%s\t%s\t%s\n' "${label}" "${assembly}" "${read1}" "${read2}" >> "${READ_MANIFEST}"
done < "${MANIFEST}"

while IFS=$'\t' read -r label assembly read1 read2; do
    [[ "${label}" == "label" ]] && continue
    sample_dir="${OUT_ROOT}/read_mapping/per_sample/${label}"
    mkdir -p "${sample_dir}"
    bam="${STAGING}/${label}.sorted.bam"
    echo "Mapping reads for ${label}" | tee -a "${OUT_ROOT}/logs/read_mapping.log"
    "${MINIMAP2}" -ax sr -t 16 "${assembly}" "${read1}" "${read2}" 2> "${sample_dir}/minimap2.log" \
        | "${SAMTOOLS}" sort -@ 4 -m 1G -o "${bam}" - 2> "${sample_dir}/samtools_sort.log"
    "${SAMTOOLS}" index -@ 4 "${bam}"
    "${SAMTOOLS}" flagstat -@ 4 -O tsv "${bam}" > "${sample_dir}/flagstat.tsv"
    "${SAMTOOLS}" stats -@ 4 "${bam}" > "${sample_dir}/samtools_stats.txt"
    "${SAMTOOLS}" coverage "${bam}" > "${sample_dir}/coverage_by_contig.tsv"
    "${SAMTOOLS}" depth -aa "${bam}" \
        | awk 'BEGIN {n=0; covered=0; sum=0} {n++; sum+=$3; if ($3>0) covered++} END {printf "positions\tcovered_positions\tcoverage_breadth_percent\tmean_depth\n%d\t%d\t%.6f\t%.6f\n", n, covered, (n ? 100*covered/n : 0), (n ? sum/n : 0)}' \
        > "${sample_dir}/depth_summary.tsv"
    rm -f -- "${bam}" "${bam}.bai"
done < "${READ_MANIFEST}"

env \
    PATH="${BUSCO_ENV}/bin:/usr/local/bin:/usr/bin:/bin" \
    "${BUSCO_ENV}/bin/python3.12" "${BUSCO_ENV}/bin/busco" \
    --in "${STAGING}" \
    --out terreus_genome_busco \
    --out_path "${OUT_ROOT}/busco_results" \
    --mode genome \
    --lineage_dataset eurotiales_odb12.2 \
    --datasets_version odb12.2 \
    --download_path "${OUT_ROOT}/busco_downloads" \
    --cpu 20 \
    --miniprot \
    --scaffold_composition \
    --tar \
    --opt-out-run-stats \
    2>&1 | tee "${OUT_ROOT}/logs/busco.log"

python3 "${OUT_ROOT}/scripts/summarize_qc.py" "${OUT_ROOT}"
date --iso-8601=seconds > "${OUT_ROOT}/QC_COMPLETE.txt"
