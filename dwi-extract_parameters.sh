#!/usr/bin/env bash

########################################
# CONFIGURATION
########################################

jsonfile=$1

# Bin width (s/mm^2) used to group b-values that should represent the same
# shell but differ slightly due to scanner rounding (e.g. 999.98 vs 1000).
BVAL_BIN_WIDTH=20

# Colors
NC='\033[0m'
BLUE='\033[0;34m'
GREEN='\033[0;32m'
RED='\033[0;31m'

Usage() {
  echo "Usage: $0 <path_to_spec.json>"
  echo "  check and store scan parameters"
  exit 1
}

log() {
    local color="$1"
    shift
    echo -e "${color}$*${NC}"
}

if [[ -z "$jsonfile" || ! -f "$jsonfile" ]]; then
  log "$RED" "Error: spec.json file not provided or does not exist."
  Usage
fi

for key in $(jq -r 'keys[]' "${jsonfile}"); do
  value=$(jq -r --arg k "$key" '.[$k]' "${jsonfile}")
  declare "$key"="$value"
done

for var in bidsdir outputdir ; do
  if [[ -z "${!var}" ]]; then
    log "$RED" "Error: Variable '$var' is not set or is empty."
    exit 1
  fi
done  

SUMMARY_FILE="${outputdir}/params/parameter_value_frequencies.json"
REPORT_FILE="${outputdir}/params/parameter_report.html"
mkdir -p "${outputdir}/params"

########################################
# PARAMETER LISTS
########################################

# NOTE: AcquisitionMatrixPE / ReconMatrixPE added here so dwi & fmap get them
# too (previously only T1w had them). InplaneResolution is NOT read from the
# json (it's computed from the NIfTI header, see get_inplane_resolution) but
# is listed here in comments for reference; it's appended separately below.
DWI_FMAP_PARAMS=(
  "CoilString"
  "MagneticFieldStrength"
  "Manufacturer"
  "ManufacturersModelName"
  "SoftwareVersions"
  "SliceThickness"
  "SpacingBetweenSlices"
  "PhaseEncodingDirection"
  "RepetitionTime"
  "InPlanePhaseEncodingDirectionDICOM"
  "EchoTime"
  "FlipAngle"
  "MultibandAccelerationFactor"
  "ParallelAcquisitionTechnique"
  "ParallelReductionFactorInPlane"
  "PixelBandwidth"
  "EffectiveEchoSpacing"
  "TotalReadoutTime"
  "PulseSequenceName"
  "ProtocolName"
  "AcquisitionMatrixPE"
  "ReconMatrixPE"
)

# FMAP_EXTRA_PARAMS=(
#   "EchoTime1"
#   "EchoTime2"
#   "Units"
# )

T1W_PARAMS=(
  "MagneticFieldStrength"
  "Manufacturer"
  "ManufacturersModelName"
  "SoftwareVersions"
  "SliceThickness"
  "EchoTime"
  "RepetitionTime"
  "InversionTime"
  "FlipAngle"
  "CoilString"
  "AcquisitionMatrixPE"
  "ReconMatrixPE"
  "InPlanePhaseEncodingDirectionDICOM"
  "PixelBandwidth"
  "ParallelReductionFactorInPlane"
  "PercentPhaseFOV"
  "ScanningSequence"
  "SequenceVariant"
  "PulseSequenceName"
  "ProtocolName"
)

########################################
# HELPERS
########################################

get_field_or_na() {
  local file="$1"
  local field="$2"

  if [[ ! -f "$file" ]]; then
    echo "NA"
    return
  fi

  local value
  value=$(jq -r --arg f "$field" '
    if has($f) and .[$f] != null then .[$f] else "NA" end
  ' "$file" 2>/dev/null)

  if [[ -z "$value" || "$value" == "null" ]]; then
    echo "NA"
  else
    echo "$value"
  fi
}

# Given a *.json sidecar path, find the matching *.nii.gz or *.nii image.
find_nii_for_json() {
  local json="$1"
  local base="${json%.json}"
  if [[ -f "${base}.nii.gz" ]]; then
    echo "${base}.nii.gz"
  elif [[ -f "${base}.nii" ]]; then
    echo "${base}.nii"
  else
    echo ""
  fi
}

# In-plane resolution isn't reliably present in the BIDS json sidecar, so
# read it straight from the NIfTI-1 header (pixdim[1] x pixdim[2], the voxel
# size along the two in-plane axes) using only python3's stdlib (struct +
# gzip), so this works whether the image is .nii or .nii.gz and without
# depending on nibabel/FSL being installed.
get_inplane_resolution() {
  local nii_file="$1"

  if [[ -z "$nii_file" || ! -f "$nii_file" ]]; then
    echo "NA"
    return
  fi

  python3 - "$nii_file" 2>/dev/null << 'PYEOF'
import sys, struct, gzip

path = sys.argv[1]

def read_header(p):
    opener = gzip.open if p.endswith(".gz") else open
    with opener(p, "rb") as f:
        return f.read(348)

def get_pixdim(data):
    # NIfTI-1 header: sizeof_hdr (int32) must read as 348. Try both
    # endiannesses since files can be written either way.
    for endian in ("<", ">"):
        if len(data) < 348:
            continue
        sizeof_hdr = struct.unpack(endian + "i", data[0:4])[0]
        if sizeof_hdr == 348:
            # pixdim is float32[8] starting at byte offset 76
            return struct.unpack(endian + "8f", data[76:76 + 32])
    return None

try:
    data = read_header(path)
    pixdim = get_pixdim(data)
    if pixdim is None:
        print("NA")
    else:
        x, y = pixdim[1], pixdim[2]
        print(f"{x:.4f}x{y:.4f}mm")
except Exception:
    print("NA")
PYEOF
}

# Compact per-run summary of exact b-values present, e.g. "b0=7, 1000=64, 2000=64"
get_bval_summary() {
  local dwi_json="$1"
  local bval_file="${dwi_json%.json}.bval"

  if [[ ! -f "$bval_file" ]]; then
    echo "NA"
    return
  fi

  tr -s ' \t' '\n' < "$bval_file" | awk -v w="$BVAL_BIN_WIDTH" '
    NF {
      v = $1 + 0
      b = int((v + w/2) / w) * w
      count[b]++
    }
    END {
      n = 0
      for (b in count) order[n++] = b
      for (i = 0; i < n; i++)
        for (j = i + 1; j < n; j++)
          if (order[j] + 0 < order[i] + 0) { t = order[i]; order[i] = order[j]; order[j] = t }
      out = ""
      for (i = 0; i < n; i++) {
        label = (order[i] == 0) ? "b0" : order[i]
        if (out != "") out = out ", "
        out = out label "=" count[order[i]]
      }
      print out
    }
  '
}

# Logs every individual b-value (rounded) as its own row so the sample-wide
# frequency summary/report shows exact counts, including b0, across all
# subjects AND all runs.
log_bvals() {
  local dwi_json="$1"
  local bval_file="${dwi_json%.json}.bval"
  [[ -f "$bval_file" ]] || return

  tr -s ' \t' '\n' < "$bval_file" | awk -v w="$BVAL_BIN_WIDTH" '
    NF {
      v = $1 + 0
      b = int((v + w/2) / w) * w
      printf "dwi\tbvalue\t%d\n", b
    }
  ' >> "${outputdir}/params/.param_values_tmp.tsv"
}

escape_json_string() {
  local s="$1"
  s="${s//\\/\\\\}"
  s="${s//\"/\\\"}"
  echo "$s"
}

log_param_value() {
  local modality="$1"
  local param="$2"
  local value="$3"

  value="${value//$'\n'/ }"
  echo -e "${modality}\t${param}\t${value}" >> "${outputdir}/params/.param_values_tmp.tsv"
}

# Build a JSON array of per-run parameter objects for one modality.
# Each element is tagged with the source filename so multiple runs
# (e.g. dir-AP / dir-PA DWI, or multiple fmap runs) are kept distinct
# rather than one silently overwriting/hiding the other.
#
# Args: modality_label  json_file_1 [json_file_2 ...]
# Reads the param list from the global array named in $PARAM_LIST_NAME.
build_modality_array() {
  local modality="$1"; shift
  local -n params_ref="$PARAM_LIST_NAME"
  local files=("$@")

  if [[ ${#files[@]} -eq 0 ]]; then
    echo "\"NA\""
    return
  fi

  local arr="[" first=1
  for f in "${files[@]}"; do
    local -A p=()

    for prm in "${params_ref[@]}"; do
      val=$(get_field_or_na "$f" "$prm")
      p["$prm"]="$val"
      log_param_value "$modality" "$prm" "$val"
    done

    # InplaneResolution: not in the json, read from the paired NIfTI header.
    nii_file="$(find_nii_for_json "$f")"
    inplane="$(get_inplane_resolution "$nii_file")"
    p["InplaneResolution"]="$inplane"
    log_param_value "$modality" "InplaneResolution" "$inplane"

    # DWI-specific: b-value shell breakdown.
    if [[ "$modality" == "dwi" ]]; then
      bval_summary=$(get_bval_summary "$f")
      p["DiffusionScheme"]="$bval_summary"
      log_param_value "dwi" "DiffusionScheme" "$bval_summary"
      log_bvals "$f"
    fi

    entry="{\"file\": \"$(basename "$f")\""
    for k in "${!p[@]}"; do
      v="$(escape_json_string "${p[$k]}")"
      entry+=", \"$k\": \"$v\""
    done
    entry+="}"

    [[ $first -eq 0 ]] && arr+=", "
    arr+="$entry"
    first=0
  done
  arr+="]"
  echo "$arr"
}

########################################
# MAIN LOOP: iterate subjects & sessions
########################################

> "${outputdir}/params/.param_values_tmp.tsv"

find "$bidsdir" -maxdepth 1 -type d -name "sub-*" | sort | while read -r SUBDIR; do
  SUBID=$(basename "$SUBDIR")

  SESSIONS=()
  while IFS= read -r ses; do
    SESSIONS+=("$ses")
  done < <(find "$SUBDIR" -maxdepth 1 -type d -name "ses-*" | sort)

  if [[ ${#SESSIONS[@]} -eq 0 ]]; then
    SESSIONS=("$SUBDIR")
  fi

  for SESDIR in "${SESSIONS[@]}"; do
    SESSION_LABEL=""
    if [[ "$SESDIR" != "$SUBDIR" ]]; then
      SESSION_LABEL=$(basename "$SESDIR")
      log "$BLUE" "Processing $SUBID $SESSION_LABEL"
    else
      log "$BLUE" "Processing $SUBID"
    fi

    DWI_DIR="$SESDIR/dwi"
    FMAP_DIR="$SESDIR/fmap"
    T1W_DIR="$SESDIR/anat"

    # Collect ALL matching sidecars, not just the first one alphabetically.
    # This is what lets multiple DWI runs with opposite phase-encoding
    # directions (e.g. dir-AP / dir-PA, or run-1 / run-2) each get their
    # own entry instead of one silently shadowing the other.
    DWI_JSONS=()
    while IFS= read -r f; do DWI_JSONS+=("$f"); done < <(find "$DWI_DIR" -maxdepth 1 -type f -name "*_dwi.json" 2>/dev/null | sort)

    FMAP_JSONS=()
    if [[ -d "$FMAP_DIR" ]]; then
      while IFS= read -r f; do FMAP_JSONS+=("$f"); done < <(find "$FMAP_DIR" -maxdepth 1 -type f -name "*_epi.json" 2>/dev/null | sort)
    fi

    T1W_JSONS=()
    while IFS= read -r f; do T1W_JSONS+=("$f"); done < <(find "$T1W_DIR" -maxdepth 1 -type f -name "*_T1w.json" 2>/dev/null | sort)

    if [[ ${#DWI_JSONS[@]} -eq 0 ]]; then
      log "$RED" "  No DWI data found for ${SUBID}${SESSION_LABEL:+ $SESSION_LABEL}."
    elif [[ ${#DWI_JSONS[@]} -gt 1 ]]; then
      log "$BLUE" "  Found ${#DWI_JSONS[@]} DWI runs for ${SUBID}${SESSION_LABEL:+ $SESSION_LABEL} (each logged separately)."
    fi

    if [[ ${#FMAP_JSONS[@]} -eq 0 ]]; then
      log "$RED" "  No fmap data for ${SUBID}${SESSION_LABEL:+ $SESSION_LABEL}. Setting fmap to \"NA\" and skipping fmap parameters."
    elif [[ ${#FMAP_JSONS[@]} -gt 1 ]]; then
      log "$BLUE" "  Found ${#FMAP_JSONS[@]} fmap runs for ${SUBID}${SESSION_LABEL:+ $SESSION_LABEL} (each logged separately)."
    fi

    PARAM_LIST_NAME="DWI_FMAP_PARAMS"
    dwi_json="$(build_modality_array "dwi" "${DWI_JSONS[@]}")"

    PARAM_LIST_NAME="DWI_FMAP_PARAMS"
    fmap_json="$(build_modality_array "fmap" "${FMAP_JSONS[@]}")"

    PARAM_LIST_NAME="T1W_PARAMS"
    t1w_json="$(build_modality_array "T1w" "${T1W_JSONS[@]}")"

    OUTPUT_JSON="{"
    OUTPUT_JSON+="\"subject\": \"${SUBID}\""
    if [[ -n "$SESSION_LABEL" ]]; then
      OUTPUT_JSON+=", \"session\": \"${SESSION_LABEL}\""
    fi
    OUTPUT_JSON+=", \"dwi\": ${dwi_json}"
    OUTPUT_JSON+=", \"fmap\": ${fmap_json}"
    OUTPUT_JSON+=", \"T1w\": ${t1w_json}"
    OUTPUT_JSON+="}"

    if [[ -n "$SESSION_LABEL" ]]; then
      OUTFILE="${outputdir}/params/${SUBID}_${SESSION_LABEL}_parameters.json"
    else
      OUTFILE="${outputdir}/params/${SUBID}_parameters.json"
    fi

    echo "$OUTPUT_JSON" | jq '.' > "$OUTFILE"

    log "$GREEN" "  Written: $OUTFILE"
  done
done

########################################
# BUILD FREQUENCY SUMMARY WITH jq
########################################

if [[ ! -s "${outputdir}/params/.param_values_tmp.tsv" ]]; then
  log "$RED" "No parameter values were logged; frequency summary will be empty."
  echo '{}' > "$SUMMARY_FILE"
else
  # Operating directly on "." (rather than a bound variable) lets jq auto-vivify
  # nested objects, and jq treats `null + 1` as `1`, so a plain reduce with +=
  # does everything the earlier has()-check version was trying to do, without
  # the "Invalid path expression" error that comes from assigning into a variable.
  awk -F '\t' 'NF==3 {print}' "${outputdir}/params/.param_values_tmp.tsv" \
    | jq -Rn '
      reduce ( inputs | split("\t") | {modality: .[0], param: .[1], value: .[2]} ) as $row (
        {};
        .[$row.modality][$row.param][$row.value] += 1
      )
    ' > "$SUMMARY_FILE"
fi

########################################
# BUILD HTML OVERVIEW REPORT (digestible, self-contained, no internet needed)
########################################

JQ_REPORT_SCRIPT="$(mktemp)"
cat > "$JQ_REPORT_SCRIPT" << 'JQEOF'
def esc: tostring | @html;

# NOTE: parameters here are declared with "$" (e.g. $pname, not pname).
# Plain (non-$) jq function params are lazy filters that get RE-EVALUATED
# against whatever "." happens to be at the point of use -- so pname (a
# filter like ".key") silently picked up the wrong .key once evaluated
# inside a nested map() with a different current input. "$"-params are
# bound to a concrete value at call time, which avoids that trap.
def bvalue_label($pname; $key):
  if $pname == "bvalue" then
    (if $key == "0" then "b0" else ($key + " s/mm\u00b2") end)
  else
    $key
  end;

def bar_row($pname; $key; $count; $maxcount):
  (if $maxcount == 0 then 0 else (($count / $maxcount) * 100) end) as $pct
  | "<div class=\"bar-row\"><span class=\"bar-label\">" + (bvalue_label($pname; $key) | esc)
    + "</span><span class=\"bar-track\"><span class=\"bar-fill\" style=\"width:" + ($pct|tostring) + "%\"></span></span>"
    + "<span class=\"bar-count\">" + ($count|tostring) + "</span></div>";

def param_section($pname; $pdata):
  ([$pdata[]] | max) as $maxcount
  | ($pdata | length) as $nvals
  | "<div class=\"param " + (if $nvals > 1 then "varies" else "consistent" end) + "\">"
  + "<h4>" + ($pname|esc) + " "
  + (if $nvals > 1 then "<span class=\"flag warn\">varies &ndash; " + ($nvals|tostring) + " values</span>"
     else "<span class=\"flag ok\">consistent</span>" end)
  + "</h4>"
  + ( $pdata | to_entries | sort_by(-.value) | map(bar_row($pname; .key; .value; $maxcount)) | join("") )
  + "</div>";

def modality_section($mname; $mdata):
  "<section><h2>" + ($mname|esc) + "</h2>"
  + ( $mdata | to_entries | sort_by(if .key == "bvalue" then 0 else 1 end)
      | map(param_section(.key; .value)) | join("") )
  + "</section>";

"<!DOCTYPE html><html><head><meta charset=\"utf-8\"><title>Scan parameter report</title><style>"
+ "body{font-family:-apple-system,Helvetica,Arial,sans-serif;max-width:900px;margin:2rem auto;padding:0 1rem;color:#1a1a1a;background:#fafafa;}"
+ "h1{font-size:1.4rem;} h2{border-bottom:2px solid #ddd;padding-bottom:.3rem;margin-top:2.5rem;} h4{margin:1rem 0 .3rem;font-size:.95rem;}"
+ ".param{padding:.5rem .75rem;margin-bottom:.4rem;border-left:4px solid #ccc;background:#fff;border-radius:4px;}"
+ ".param.varies{border-left-color:#d9534f;background:#fff7f6;}"
+ ".param.consistent{border-left-color:#5cb85c;}"
+ ".flag{font-size:.7rem;font-weight:600;padding:.1rem .4rem;border-radius:3px;margin-left:.3rem;}"
+ ".flag.warn{background:#f8d7da;color:#842029;} .flag.ok{background:#d1e7dd;color:#0f5132;}"
+ ".bar-row{display:flex;align-items:center;gap:.5rem;font-size:.82rem;margin:.15rem 0;}"
+ ".bar-label{flex:0 0 140px;text-align:right;color:#333;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;}"
+ ".bar-track{flex:1;background:#eee;border-radius:3px;height:14px;position:relative;}"
+ ".bar-fill{display:block;height:100%;background:#4a90d9;border-radius:3px;min-width:2px;}"
+ ".param.varies .bar-fill{background:#d9534f;}"
+ ".bar-count{flex:0 0 30px;color:#666;}"
+ "p.note{color:#666;font-size:.85rem;}"
+ "</style></head><body>"
+ "<h1>Scan parameter overview</h1>"
+ "<p class=\"note\">Green = consistent across the sample. Red = varies &ndash; check whether this reflects a real protocol difference or a header/conversion issue.</p>"
+ ( to_entries | map(modality_section(.key; .value)) | join("") )
+ "</body></html>"
JQEOF

jq -r -f "$JQ_REPORT_SCRIPT" "$SUMMARY_FILE" > "$REPORT_FILE"
rm -f "$JQ_REPORT_SCRIPT"

########################################
# DISPLAY SUMMARY OF DIFFERENCES
########################################

echo
log "$BLUE" "Parameter value frequencies (JSON) summarized in:"
echo "  $SUMMARY_FILE"
jq '.' "$SUMMARY_FILE"

echo
log "$BLUE" "Human-readable summary of frequencies:"
jq -r '
  to_entries[] as $m |
  "Modality: \($m.key)" ,
  (
    $m.value
    | to_entries[]
    | "  Parameter: \(.key)",
      (
        .value
        | to_entries[]
        | "    \(.key): \(.value) occurrences"
      )
  ),
  ""
' "$SUMMARY_FILE"

echo
log "$GREEN" "If a parameter has more than one value for a modality, it indicates differences in the dataset."
echo
log "$BLUE" "Visual overview (bar charts, flags parameters that vary):"
echo "  $REPORT_FILE"
log "$GREEN" "Open it in a browser for a quick visual check across the sample."