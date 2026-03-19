#!/usr/bin/env bash
set -euo pipefail

# ---- CONFIG ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GEF_SCRIPT="$SCRIPT_DIR/GEFbashscript_nolmd.sh"

Z=92
A=239
FACTOR=1

SPIN_MIN=0.5
SPIN_MAX=17.5
SPIN_STEP=1

E_MIN=4.6
E_MAX=24
E_STEP=0.2

echo "Starting GEF grid run"
echo "Z=$Z A=$A  Spin $SPIN_MIN→$SPIN_MAX  E* $E_MIN→$E_MAX step $E_STEP"
echo

# -------- LOOP ON SPIN --------
awk -v Smin="$SPIN_MIN" -v Smax="$SPIN_MAX" -v Sstep="$SPIN_STEP" '
BEGIN {
  eps = Sstep/10.0;
  for (S = Smin; S <= Smax + eps; S += Sstep) {
    printf "%.1f\n", S
  }
}' | while read -r spin; do

  echo "===== Spin = $spin ====="

  # -------- LOOP ON E* --------
  awk -v Emin="$E_MIN" -v Emax="$E_MAX" -v Estep="$E_STEP" '
  BEGIN {
    eps = Estep/10.0;
    for (E = Emin; E <= Emax + eps; E += Estep) {
      printf "%.2f\n", E
    }
  }' | while read -r EEXC; do

    echo "Running E* = $EEXC  Spin = $spin"

    "$GEF_SCRIPT" "$Z" "$A" "$EEXC" "$spin" "$FACTOR"

  done

done

echo
echo "All runs finished CHECK"