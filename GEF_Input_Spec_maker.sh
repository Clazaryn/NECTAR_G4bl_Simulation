#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./GEF_Input_Spec_maker.sh <GEF_DIR> <EEXC> <EEXC_NEXT> <STEP>
#
# Produces:
#   <GEF_DIR>/in/StepSpectrum_<EEXC>_<EEXC_NEXT>
# with 2 columns:
#   E (from EEXC to EEXC_NEXT inclusive, step STEP)   0.1

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <GEF_DIR> <EEXC> <EEXC_NEXT> <STEP>" >&2
  exit 1
fi

GEF_DIR="$1"
EEXC="$2"
EEXC_NEXT="$3"
STEP="$4"

IN_DIR="$GEF_DIR/in"
mkdir -p "$IN_DIR"

SPEC_BASENAME="StepSpectrum_${EEXC}_${EEXC_NEXT}_${STEP}.in"
SPEC_PATH="$IN_DIR/$SPEC_BASENAME"

# If already exists and non-empty, do nothing
if [[ -s "$SPEC_PATH" ]]; then
  echo "$SPEC_PATH"
  exit 0
fi

# Generate spectrum
# Note: use LC_NUMERIC=C to ensure dot decimal separator
LC_NUMERIC=C awk -v s="$EEXC" -v e="$EEXC_NEXT" -v step="$STEP" '
BEGIN {
  s += 0; e += 0; step += 0;
  if (step <= 0) { exit 2 }
  # include end point safely with epsilon
  eps = step/10.0;
  for (x = s; x <= e + eps; x += step) {
    printf "%.2f 0.1\n", x
  }
}
' > "$SPEC_PATH"

echo "$SPEC_PATH"