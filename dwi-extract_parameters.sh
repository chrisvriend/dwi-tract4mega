#!/usr/bin/env bash

########################################
# CONFIGURATION
########################################

jsonfile=$1

# Bin width (s/mm^2) used to group b-values that should represent the same
# shell but differ slightly due to scanner rounding (e.g. 999.98 vs 1000).
BVAL_BIN_WIDTH=50

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

SUMMARY_FILE="${outputdir}/tosend/parameter_value_frequencies.json"
REPORT_FILE="${outputdir}/tosend/parameter_report.html"
mkdir -p "${outputdir}/tosend"

########################################
# PARAMETER LISTS
########################################

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
)

FMAP_EXTRA_PARAMS=(
  "EchoTime1"
  "EchoTime2"
  "Units"
)

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

# Compact per-subject summary of exact b-values present, e.g. "b0=7, 1000=64, 2000=64"
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
# frequency summary/report shows exact counts, including b0, across all subjects.
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
  ' >> "${outputdir}/tosend/.param_values_tmp.tsv"
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
  echo -e "${modality}\t${param}\t${value}" >> "${outputdir}/tosend/.param_values_tmp.tsv"
}

########################################
# MAIN LOOP: iterate subjects & sessions
########################################

> "${outputdir}/tosend/.param_values_tmp.tsv"

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

    DWI_JSON="$(find "$DWI_DIR" -maxdepth 1 -type f -name "*_dwi.json" 2>/dev/null | sort | head -n 1)"

    FMAP_JSON=""
    if [[ -d "$FMAP_DIR" ]]; then
      FMAP_JSON="$(find "$FMAP_DIR" -maxdepth 1 -type f -name "*_epi.json" 2>/dev/null | sort | head -n 1)"
    fi

    T1W_JSON="$(find "$T1W_DIR" -maxdepth 1 -type f -name "*_T1w.json" 2>/dev/null | sort | head -n 1)"

    # DWI params
    declare -A dwi_params
    for p in "${DWI_FMAP_PARAMS[@]}"; do
      val=$(get_field_or_na "${DWI_JSON:-/dev/null}" "$p")
      dwi_params["$p"]="$val"
      log_param_value "dwi" "$p" "$val"
    done

    # b-value breakdown (per-subject compact string + individual rows for aggregate counting)
    bval_summary=$(get_bval_summary "${DWI_JSON:-}")
    dwi_params["DiffusionScheme"]="$bval_summary"
    log_param_value "dwi" "DiffusionScheme" "$bval_summary"
    log_bvals "${DWI_JSON:-}"

    # fmap params
    declare -A fmap_params
    fmap_available=true
    if [[ ! -d "$FMAP_DIR" || -z "$FMAP_JSON" ]]; then
      fmap_available=false
      if [[ -n "$SESSION_LABEL" ]]; then
        log "$RED" "  No fmap data for ${SUBID} ${SESSION_LABEL}. Setting fmap to \"NA\" and skipping fmap parameters."
      else
        log "$RED" "  No fmap data for ${SUBID}. Setting fmap to \"NA\" and skipping fmap parameters."
      fi
    fi

    if $fmap_available; then
      for p in "${DWI_FMAP_PARAMS[@]}"; do
        val=$(get_field_or_na "${FMAP_JSON:-/dev/null}" "$p")
        fmap_params["$p"]="$val"
        log_param_value "fmap" "$p" "$val"
      done
      for p in "${FMAP_EXTRA_PARAMS[@]}"; do
        val=$(get_field_or_na "${FMAP_JSON:-/dev/null}" "$p")
        fmap_params["$p"]="$val"
        log_param_value "fmap" "$p" "$val"
      done
    fi

    # T1w params
    declare -A t1w_params
    for p in "${T1W_PARAMS[@]}"; do
      val=$(get_field_or_na "${T1W_JSON:-/dev/null}" "$p")
      t1w_params["$p"]="$val"
      log_param_value "T1w" "$p" "$val"
    done

    # Build dwi JSON
    dwi_json="{"; first=1
    for k in "${!dwi_params[@]}"; do
      v="$(escape_json_string "${dwi_params[$k]}")"
      [[ $first -eq 0 ]] && dwi_json+=", "
      dwi_json+="\"$k\": \"$v\""
      first=0
    done
    dwi_json+="}"

    # Build fmap JSON
    if $fmap_available; then
      fmap_json="{"; first=1
      for k in "${!fmap_params[@]}"; do
        v="$(escape_json_string "${fmap_params[$k]}")"
        [[ $first -eq 0 ]] && fmap_json+=", "
        fmap_json+="\"$k\": \"$v\""
        first=0
      done
      fmap_json+="}"
    else
      fmap_json="\"NA\""
    fi

    # Build T1w JSON
    t1w_json="{"; first=1
    for k in "${!t1w_params[@]}"; do
      v="$(escape_json_string "${t1w_params[$k]}")"
      [[ $first -eq 0 ]] && t1w_json+=", "
      t1w_json+="\"$k\": \"$v\""
      first=0
    done
    t1w_json+="}"

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
      OUTFILE="${outputdir}/tosend/${SUBID}_${SESSION_LABEL}_parameters.json"
    else
      OUTFILE="${outputdir}/tosend/${SUBID}_parameters.json"
    fi

    echo "$OUTPUT_JSON" | jq '.' > "$OUTFILE"

    log "$GREEN" "  Written: $OUTFILE"

    unset dwi_params fmap_params t1w_params
  done
done

########################################
# BUILD FREQUENCY SUMMARY WITH jq
########################################

if [[ ! -s "${outputdir}/tosend/.param_values_tmp.tsv" ]]; then
  log "$RED" "No parameter values were logged; frequency summary will be empty."
  echo '{}' > "$SUMMARY_FILE"
else
  # Operating directly on "." (rather than a bound variable) lets jq auto-vivify
  # nested objects, and jq treats `null + 1` as `1`, so a plain reduce with +=
  # does everything the earlier has()-check version was trying to do, without
  # the "Invalid path expression" error that comes from assigning into a variable.
  awk -F '\t' 'NF==3 {print}' "${outputdir}/tosend/.param_values_tmp.tsv" \
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