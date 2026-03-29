#!/bin/bash -l
#SBATCH -J irredeemable_cleanup
#SBATCH -p cloud
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o sbatch-irredeemable-cleanup-%j.out
#SBATCH -e sbatch-irredeemable-cleanup-%j.err

set -euo pipefail

repo_root=${SNHEART_ROOT:-/path/to/SNHEART}
analysis_dir="$repo_root/MRI_data/analysis"
cleanup_ids=(
  546 564 569 575 577 580 581 596 615 616 617
  619 620 621 628 631 632 633 634 635 637
)
cleanup_joined=$(IFS=,; echo "${cleanup_ids[*]}")

mapfile -t thr_files < <(ls -1 "$analysis_dir"/SN_nmdata_irredeemable_*_thr10.txt 2>/dev/null | sort)
if [[ ${#thr_files[@]} -eq 0 ]]; then
  echo "[ERROR] No existing SN_nmdata_irredeemable_*_thr10.txt table found in $analysis_dir" >&2
  exit 1
fi
current_thr="${thr_files[-1]}"

echo "Current Irredeemable table: $current_thr"
existing_cleanup_rows=$(python - <<'PY' "$current_thr" "$cleanup_joined"
import csv, sys
cur = sys.argv[1]
ids = set(sys.argv[2].split(','))
count = 0
with open(cur) as fin:
    reader = csv.DictReader(fin)
    for row in reader:
        if row['sub'] in ids:
            count += 1
print(count)
PY
)
echo "Cleanup rows already present in current table: $existing_cleanup_rows"

export MATLAB_ROOT=${MATLAB_ROOT:-"$repo_root/Software"}
MATLAB_BIN="$MATLAB_ROOT/bin/matlab"
if [[ ! -x $MATLAB_BIN ]]; then
  echo "[ERROR] MATLAB binary not found at $MATLAB_BIN" >&2
  exit 1
fi

log_dir="$analysis_dir/logs"
mkdir -p "$log_dir"
job_tag=${SLURM_JOB_ID:-$$}
export IRREDEEMABLE_LOG="$log_dir/irredeemable_cleanup_${job_tag}.log"
export IRREDEEMABLE_SUBJECTS="${cleanup_ids[*]}"

"$MATLAB_BIN" -nodisplay -nosplash -nodesktop \
  -r "addpath('$repo_root/Scripts'); run('Scripts/irredeemable_black_ant_V1.m'); exit;"

echo "MATLAB run complete. Appending cleanup rows to $current_thr"

mapfile -t new_thr_files < <(ls -1 "$analysis_dir"/SN_nmdata_irredeemable_*_thr10.txt 2>/dev/null | sort)
new_thr="${new_thr_files[-1]}"
if [[ $new_thr == "$current_thr" ]]; then
  echo "[ERROR] No new Irredeemable table detected after MATLAB run" >&2
  exit 1
fi

tmp_filtered="$current_thr.filtered"
removed_rows=$(python - <<'PY' "$current_thr" "$tmp_filtered" "$cleanup_joined"
import csv, sys
cur, tmp, ids = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
cleanup = set(ids)
removed = 0
with open(cur) as fin:
    reader = csv.reader(fin)
    header = next(reader)
    rows = []
    for row in reader:
        if row and row[0] in cleanup:
            removed += 1
        else:
            rows.append(row)
with open(tmp, 'w', newline='') as fout:
    writer = csv.writer(fout)
    writer.writerow(header)
    writer.writerows(rows)
print(removed)
PY
)
mv "$tmp_filtered" "$current_thr"
echo "Removed $removed_rows existing cleanup rows from $current_thr"

rows_to_add=$(($(wc -l < "$new_thr") - 1))
if (( rows_to_add <= 0 )); then
  echo "[ERROR] New table $new_thr has no data rows to append" >&2
  exit 1
fi

tail -n +2 "$new_thr" >> "$current_thr"
echo "Appended $rows_to_add cleanup rows from $new_thr to $current_thr"

echo "Done. Review $current_thr and consider archiving $new_thr if desired."
