#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_gef.sh 92 238 9.0 10
# Produces:
#   GEFResults_Z92_A238_E9.0_factor_10.lmd   (in the GEF dir by default)

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <Z> <A> <Eexc_MeV> <EexcNext_MeV> <enhancement_factor>"
  exit 1
fi

Z="$1"
A="$2"
EEXC="$3"
EEXC_NEXT="$4"
FACTOR="$5"

# --- Configure where GEF lives ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GEF_DIR="../GEF"   # TODO: set this
GEF_BIN="./GEF"               # or ./GEF.exe etc.

# ---------------------------------
# Ensure EEXC spectrum file exists
# ---------------------------------
STEP="0.01"
SPEC_BASENAME="StepSpectrum_${EEXC}_${EEXC_NEXT}_${STEP}.in"
SPEC_PATH="$GEF_DIR/in/$SPEC_BASENAME"

if [[ ! -s "$SPEC_PATH" ]]; then
  echo "Spectrum file missing -> generating: $SPEC_PATH"
  "$SCRIPT_DIR/GEF_Input_Spec_maker.sh" "$GEF_DIR" "$EEXC" "$EEXC_NEXT" "$STEP" >/dev/null
fi

# This is what we pass to GEF (relative to GEF_DIR is usually safest)
EEXC_SPEC_FILE="in/$SPEC_BASENAME"


# --- Output file name (exact format you requested) ---
EEXC_FMT=$(LC_NUMERIC=C printf "%.1f" "$EEXC") #Forces the format so it can be read by the GEF reader after...
LMD_FILE="TEST_GEFResults_RandinBin_Z${Z}_A${A}_E${EEXC_FMT}_factor_${FACTOR}.lmd"

# --- Inputs to feed to GEF ---
# IMPORTANT: do NOT put '#' comments inside the input stream.
inputs=(
  "no"            # Output of FY in ENDF format required? [y/n]
  "${Z} ${A}"     # Z and A of fissioning nucleus
  "EM"            # energy input option (as in your screenshot)
  "${EEXC_SPEC_FILE}"       # spectrum file

  ""              # ENTER
  ""              # ENTER
  ""              # ENTER
  ""              # ENTER
  ""              # ENTER

  "${FACTOR}"     # Enhancement factor
  "${LMD_FILE}"   # list-mode output filename

  ""              # ENTER (neutrons list? default no)
  ""              # ENTER (gammas list? default no)
  ""              # ENTER (Eexc contributions? default no)

  "n"            # Show display of mass distributions? [y/n]
)

echo "MAKING $GEF_DIR/$LMD_FILE IT TAKE SOME TIME, BE PATIENT !"

# First check if it the file already exist if not then it does the job if not !
OUTPUT_PATH="$GEF_DIR/out/$LMD_FILE"

if [ -s "$OUTPUT_PATH" ]; then
    echo "Skipping $LMD_FILE (already computed)"
else
    (
      cd "$GEF_DIR"
      printf "%s\n" "${inputs[@]}" | "$GEF_BIN" > /dev/null
    )
fi

# ---------------------------------
# Run ROOT macro only if tree missing
# ---------------------------------
#
#TREE_FILE="./GEF_tree/GEFResults_Z${Z}_A${A}_E${EEXC_FMT}_factor_${FACTOR}.root"
#
#if [[ -s "$TREE_FILE" ]]; then
#  echo "Skipping ROOT tree (already exists): $TREE_FILE"
#else
#  echo "Building ROOT tree: $TREE_FILE"
#  # Run from SCRIPT_DIR so relative paths inside the macro (like ./out/...) land here
#  (
#    cd "$SCRIPT_DIR"
#    root -l -b -q "./UtilityScripts/GEF_Reader.cpp+($Z, $A, $EEXC_FMT, $FACTOR);"
#  )
#
#  if [[ ! -s "$TREE_FILE" ]]; then
#    echo "WARNING: ROOT finished but tree file not found (or empty): $TREE_FILE" >&2
#    echo "         If your macro writes to a different name/location, update TREE_FILE in this script." >&2
#  fi
#fi


echo "Done. Wrote: $GEF_DIR/$LMD_FILE"

