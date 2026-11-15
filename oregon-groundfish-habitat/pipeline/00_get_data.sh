#!/usr/bin/env bash
# Fetch the NWFSC survey tables and 2-km grid from the public surveyjoin repositories.
# Usage: bash pipeline/00_get_data.sh   (writes raw/surveyjoin-data and raw/surveyjoin)
set -euo pipefail
mkdir -p raw && cd raw
fetch () {  # repo, commit, files...
  local repo=$1 commit=$2; shift 2
  [ -d "$repo" ] || git clone -q --filter=blob:none --no-checkout "https://github.com/DFO-NOAA-Pacific/$repo" "$repo"
  git -C "$repo" fetch -q origin "$commit" && git -C "$repo" checkout -q "$commit" -- "$@"
}
fetch surveyjoin-data 4f898a4cdf8f23793b0aae55a1db86b01d2f50ec nwfsc-haul.rds nwfsc-catch.rds   # v0.0.4
fetch surveyjoin      e9f4f229a30f61bd837e30ab32af608ebdde6854 data/nwfsc_grid.rda
echo "raw data ready"
