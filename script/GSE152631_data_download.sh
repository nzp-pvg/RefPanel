#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------
# Download GEO raw files for GSE152631
# English comments only.
#
# What you get:
#   1) A local folder ./GSE152631_RAW
#   2) GEO supplementary raw package(s) (usually .tar or .tar.gz)
#   3) Automatically extracted files (if tar is present)
#
# Requirements:
#   - curl (or wget)
#   - tar
# Optional but recommended:
#   - md5sum (Linux) or md5 (macOS) for checksum verification
# --------------------------------------------

GSE="GSE152631"
OUTDIR="${GSE}_RAW"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "[INFO] Downloading ${GSE} RAW supplementary files from GEO..."

# GEO stores series supplementary files under:
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSExxxxxx/suppl/
# where GSEnnn is the first 3 digits group (e.g., GSE152 -> for GSE152631)
GSE_PREFIX="${GSE:0:6}"   # e.g., GSE152
BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/${GSE_PREFIX}nnn/${GSE}/suppl/"

echo "[INFO] GEO URL: ${BASE_URL}"
echo "[INFO] Listing remote files..."

# List files available in the 'suppl' directory
# (The directory listing is plain HTML; we extract file names ending with common extensions)
curl -fsSL "${BASE_URL}" \
  | grep -Eo 'href="[^"]+"' \
  | sed -E 's/href="([^"]+)"/\1/' \
  | grep -E '\.(tar|tar\.gz|tgz|zip|gz|bz2|xz)$' \
  | sort -u \
  | tee file_list.txt

if [[ ! -s file_list.txt ]]; then
  echo "[ERROR] No downloadable files found at: ${BASE_URL}"
  echo "[HINT] Open the URL in a browser to confirm the directory contents."
  exit 1
fi

echo "[INFO] Downloading files..."
while read -r f; do
  echo "  - ${f}"
  curl -fL -O "${BASE_URL}${f}"
done < file_list.txt

echo "[INFO] Download complete."

# If MD5 files exist, verify checksums (GEO often provides .md5)
# Example: some directories include a file like "MD5SUMS" or "*.md5"
MD5_FOUND=$(ls 2>/dev/null | grep -E '(\.md5$|MD5|md5sum|MD5SUMS)' || true)
if [[ -n "${MD5_FOUND}" ]]; then
  echo "[INFO] MD5-related file(s) found:"
  echo "${MD5_FOUND}"
  echo "[INFO] If you want checksum verification, use the appropriate command:"
  echo "  - Linux: md5sum -c <md5_file>"
  echo "  - macOS: md5 -r <file>"
fi

echo "[INFO] Extracting archives when possible..."
shopt -s nullglob
for a in *.tar *.tar.gz *.tgz; do
  echo "  - extracting ${a}"
  tar -xf "${a}"
done
for z in *.zip; do
  echo "  - extracting ${z}"
  unzip -o "${z}"
done
shopt -u nullglob

echo "[DONE] ${GSE} RAW files are in: $(pwd)"
