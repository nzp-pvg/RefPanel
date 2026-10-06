#!/usr/bin/env Rscript
rm(list = ls())
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else file.path("intermediate_data", "cohort_D", "HNSCC_Zenodo_4698004")
bed_dir <- file.path(root, "bed")
files <- sort(list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE))
if (!length(files)) stop("No BED files under ", bed_dir)
bn <- basename(files)
parse <- regexec("^(cfDNA|PBL)_(control|patient)_([0-9]+)(?:_(diagnosis|BL))?\\.bed$", bn)
parts <- regmatches(bn, parse)
if (any(lengths(parts) == 0L)) stop("Unrecognized BED filename(s): ", paste(bn[lengths(parts) == 0L], collapse = ", "))
manifest <- rbindlist(Map(function(x, path) data.table(
  sample_id = sub("\\.bed$", "", basename(path)), analyte = x[2], status = x[3], subject_number = as.integer(x[4]),
  timepoint = ifelse(length(x) >= 5L && nzchar(x[5]), x[5], "baseline"), path = normalizePath(path),
  primary_cfDNA = x[2] == "cfDNA", secondary_PBL = x[2] == "PBL" && !(x[3] == "control" && x[4] == "20"),
  exclusion_reason = ifelse(x[2] == "PBL" && x[3] == "control" && x[4] == "20",
                            "Byte-for-byte identical to FaDu_MeDIP.bed", "")), parts, files))
setorder(manifest, analyte, status, subject_number)
if (manifest[primary_cfDNA == TRUE, .N] != 50L || manifest[primary_cfDNA == TRUE & status == "patient", .N] != 30L ||
    manifest[primary_cfDNA == TRUE & status == "control", .N] != 20L) stop("Unexpected cfDNA sample counts")
if (manifest[secondary_PBL == TRUE, .N] != 49L || manifest[secondary_PBL == TRUE & status == "patient", .N] != 30L ||
    manifest[secondary_PBL == TRUE & status == "control", .N] != 19L) stop("Unexpected PBL sample counts")
fwrite(manifest, file.path(root, "cohort_D_sample_manifest_v1.tsv"), sep = "\t")
message("[DONE] D manifest: cfDNA 30+20; PBL secondary 30+19; PBL_control_20 excluded")
