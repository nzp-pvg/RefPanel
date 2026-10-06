#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd "$(dirname "$0")" && pwd)"
Rscript "$script_dir/deploy_frozen_model_C_D_v1.R"
Rscript "$script_dir/score_cohort_F_frozen_model_all_v1.R"
Rscript "$script_dir/bootstrap_cohort_F_deployment_v1.R"
Rscript "$script_dir/build_transfer_benchmark_CDF_v1.R"
echo "FROZEN_DEPLOYMENT_QC_PASS"
