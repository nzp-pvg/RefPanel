#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
root="$repo_root/intermediate_data/cohort_F"
model="$repo_root/intermediate_data/frozen_model/C_model_features_v1.bed"
out="$root/derived/C_model_representation_v1"
bed_dir="${1:-$root/source/bed_full_v1}"
if [[ ! -d "$bed_dir" ]]; then
  echo "Fragment BED directory not found: $bed_dir" >&2
  echo "Download the GEO BED files listed in intermediate_data/cohort_F/source/GSE243474_MeDIP_bed_download_manifest_v1.tsv, then rerun with that directory as the first argument." >&2
  exit 2
fi
mkdir -p "$out"
for f in "$bed_dir"/*.bed.gz; do
  id=$(basename "$f" .bed.gz)
  dest="$out/$id.counts.tsv"
  [[ -s "$dest" ]] && continue
  tmp=$(mktemp)
  gzip -cd "$f" | awk 'BEGIN{OFS="\t"} NF>=3 && $2>=0 && $3>$2 && ($3-$2)>=20 && ($3-$2)<=1000 {print $1,$2,$3}' > "$tmp"
  bedtools intersect -c -a "$model" -b "$tmp" > "$dest"
  rm -f "$tmp"
  echo "$id" >&2
done
