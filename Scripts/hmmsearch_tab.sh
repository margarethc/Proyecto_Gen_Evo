#!/usr/bin/env bash
# Run HMMER (hmmsearch) with one model (*.hmm) against many protein files (*.faa / *.fasta)
# Produces one tab-delimited *_hits.txt output per protein input.
#
# Usage:
#   ./hmmsearch_tab.sh model.hmm prot1.fasta prot2.faa ...
#   ./hmmsearch_tab.sh model.hmm prot1.fasta prot2.faa /path/to/outdir/   # last arg existing directory -> outputs written there
#
# Example:
#   ./hmmsearch_tab.sh model.hmm protein_sample.fasta
#     -> produces protein_sample_hits.txt
#
# Notes:
#  - Supports .gz-compressed protein files in naming (just strips .gz when naming output).
#  - Requires hmmsearch in PATH (HMMER).
#  - Produces per-domain rows using HMMER --domtblout parsing (TSV).
set -euo pipefail

usage() {
  cat <<EOF >&2
Usage: $0 <model.hmm> <proteins.faa|.fasta[.gz]> [proteins2 ...] [out_dir]
If the last argument is an existing directory it will be used as the output directory.
EOF
  exit 2
}

if [ "$#" -lt 2 ]; then
  usage
fi

MODEL="$1"
shift

if [ ! -f "$MODEL" ]; then
  echo "Model file not found: $MODEL" >&2
  exit 1
fi

# If last arg is an existing directory, treat as outdir and remove from protein list
last_arg="${!#}"
OUTDIR=""
if [ -d "$last_arg" ]; then
  OUTDIR="${last_arg%/}"  # trim trailing slash
  # remove last arg from positional parameters
  set -- "${@:1:$(($#-1))}"
fi

# Remaining args are protein files
if [ "$#" -lt 1 ]; then
  usage
fi

# Check hmmsearch
if ! command -v hmmsearch >/dev/null 2>&1; then
  echo "hmmsearch not found in PATH. Install HMMER (http://hmmer.org/)" >&2
  exit 1
fi

# Helper: compute default output filename from protein path
default_out_from_protein() {
  local p="$1"
  local base tmpname name
  base=$(basename "$p")
  tmpname="$base"
  if [[ "$tmpname" == *.gz ]]; then
    tmpname="${tmpname%.gz}"
  fi
  name="${tmpname%.*}"
  if [ -n "$OUTDIR" ]; then
    echo "${OUTDIR}/${name}_hits.txt"
  else
    echo "${name}_hits.txt"
  fi
}

# Process each protein file
for PROTEIN in "$@"; do
  if [ ! -f "$PROTEIN" ]; then
    echo "Warning: proteins file not found, skipping: $PROTEIN" >&2
    continue
  fi

  OUT="$(default_out_from_protein "$PROTEIN")"
  TMP_DOM=$(mktemp --suffix=.domtblout)

  cleanup_file() {
    rm -f "$TMP_DOM"
  }
  trap cleanup_file EXIT

  echo "Running hmmsearch: model='$MODEL' proteins='$PROTEIN' -> out='$OUT'"

  if ! hmmsearch --domtblout "$TMP_DOM" "$MODEL" "$PROTEIN" >/dev/null 2>/dev/null; then
    echo "hmmsearch failed for $PROTEIN. Skipping." >&2
    rm -f "$TMP_DOM"
    trap - EXIT
    continue
  fi

  # Parse domtblout: skip comment lines starting with '#'
  # Output columns (tab-separated):
  # query_name, target_name, domain_iEvalue, domain_bitscore, fullseq_Evalue, fullseq_bitscore,
  # hmm_from, hmm_to, ali_from, ali_to, acc, target_description
  {
    echo -e "query_name\ttarget_name\tdomain_iEvalue\tdomain_bitscore\tfullseq_Evalue\tfullseq_bitscore\thmm_from\thmm_to\tali_from\tali_to\tacc\ttarget_description"
    grep -v '^#' "$TMP_DOM" | awk 'BEGIN{OFS="\t"} {
      desc="";
      if (NF>=23) {
        desc=$23;
        for(i=24;i<=NF;i++) desc=desc " " $i;
      }
      # HMMER domtblout columns: target name ($1), target accession ($2), tlen ($3), query name ($4), ...
      print $4, $1, $13, $14, $7, $8, $16, $17, $18, $19, $22, desc
    }'
  } > "$OUT"

  echo "Tabular results written to: $OUT"

  rm -f "$TMP_DOM"
  trap - EXIT
done
