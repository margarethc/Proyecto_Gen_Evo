#!/usr/bin/env bash
# Filter HMMER tabular hits by e-value and extract matching sequences into a new FASTA.
#
# Usage:
#   ./filter_and_extract.sh [-t THRESHOLD] [-c domain|full] [-o OUTDIR] -- hits1.txt prot1.faa [hits2.txt prot2.faa ...]
#
# Examples:
#   # default threshold 1e-10, use domain e-value (third column)
#   ./filter_and_extract.sh -- sample_hits.txt sample.faa
#
#   # explicit threshold and multiple pairs, output to ./filtered_outputs/
#   ./filter_and_extract.sh -t 1e-34 -c domain -o filtered_outputs -- s1_hits.txt s1.faa s2_hits.txt s2.faa
#
# Notes:
#  - The script expects pairs after the --: each hits file followed by its corresponding protein FASTA.
#  - The hits file must be in the tabular format produced earlier (header shown in your example).
#  - By default the script uses the "domain" e-value (column 3). Use -c full to use the full-sequence e-value (column 5).
#  - Input FASTA files can be gzipped (ending in .gz). Hits files can also be gzipped.
#  - Output FASTA is named <protein_basename>_filtered.fasta (placed in OUTDIR if provided).
set -euo pipefail

THRESHOLD="1e-10"
COLUMN="domain"   # "domain" -> use column 3, "full" -> use column 5
OUTDIR=""
PRINT_HELP=0

usage() {
  cat <<EOF >&2
Usage: $0 [-t THRESHOLD] [-c domain|full] [-o OUTDIR] -- hits1.txt prot1.faa [hits2.txt prot2.faa ...]
  -t THRESHOLD   e-value threshold (default: $THRESHOLD). Matches if evalue < THRESHOLD
  -c COLUMN      which e-value to use: domain (column 3) or full (column 5). Default: $COLUMN
  -o OUTDIR      write outputs into this directory (created if missing)
  --             end options; then give pairs: hits protein hits protein ...
EOF
  exit 2
}

while getopts ":t:c:o:h" opt; do
  case "$opt" in
    t) THRESHOLD="$OPTARG" ;;
    c) COLUMN="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    h) PRINT_HELP=1 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done
shift $((OPTIND - 1))

if [ "$PRINT_HELP" -eq 1 ]; then
  usage
fi

# require '--' to separate options from pairs (help avoid ambiguity)
if [ "$#" -lt 2 ]; then
  usage
fi

# If OUTDIR provided, create it
if [ -n "$OUTDIR" ]; then
  mkdir -p "$OUTDIR"
fi

# Determine which column index to use for e-value
case "$COLUMN" in
  domain) COL_IDX=3 ;;
  full)   COL_IDX=5 ;;
  *)
    echo "Invalid column: $COLUMN. Choose 'domain' or 'full'." >&2
    exit 2
    ;;
esac

# Ensure even number of remaining args: pairs of (hits, protein)
if [ $(( $# % 2 )) -ne 0 ]; then
  echo "Expect pairs of arguments: hits1 prot1 [hits2 prot2 ...]" >&2
  usage
fi

# helper: read possibly gzipped file via a command into stdout
read_stream_cmd() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    # prefer gzip -dc, fall back to zcat
    if command -v gzip >/dev/null 2>&1; then
      printf "gzip -dc %q" "$f"
    else
      printf "zcat %q" "$f"
    fi
  else
    printf "cat %q" "$f"
  fi
}

# Loop over pairs
while [ "$#" -gt 0 ]; do
  HITFILE="$1"; shift
  PROTEIN="$1"; shift

  if [ ! -f "$HITFILE" ]; then
    echo "Warning: hits file not found, skipping: $HITFILE" >&2
    continue
  fi
  if [ ! -f "$PROTEIN" ]; then
    echo "Warning: protein file not found, skipping: $PROTEIN" >&2
    continue
  fi

  # compute output name from protein basename
  base=$(basename "$PROTEIN")
  # strip .gz if present
  if [[ "$base" == *.gz ]]; then base="${base%.gz}"; fi
  # strip final extension
  name="${base%.*}"
  out_fasta="${name}_filtered.fasta"
  if [ -n "$OUTDIR" ]; then
    out_fasta="${OUTDIR%/}/${out_fasta}"
  fi

  echo "Processing pair:"
  echo "  hits:    $HITFILE"
  echo "  protein: $PROTEIN"
  echo "  threshold: $THRESHOLD (using column $COL_IDX - $COLUMN)"
  echo "  output:  $out_fasta"

  # build temporary id file
  idfile="$(mktemp)"
  trap 'rm -f "$idfile"' EXIT

  # Extract target names (column 2) where chosen e-value column < THRESHOLD
  # Skip header lines (which may start with "query_name" or '#')
  # Support hits file gzipped via read_stream_cmd
  eval "$(read_stream_cmd "$HITFILE")" | \
    awk -v col="$COL_IDX" -v thr="$THRESHOLD" '
      BEGIN{
        # convert threshold to numeric once (awk supports scientific notation)
        thrn = thr + 0
      }
      /^#/ { next }                             # skip HMMER comments
      NR==1 && tolower($0) ~ /query_name/ { next } # skip header line if present
      NF < col { next }                         # skip malformed lines
      {
        # numeric conversion handles scientific notation like 1e-61
        ev = ($col + 0)
        if (ev < thrn) {
          # target name is column 2
          print $2
        }
      }
    ' | sort -u > "$idfile"

  # count how many ids
  idcount=$(wc -l < "$idfile" | tr -d ' ')
  if [ "$idcount" -eq 0 ]; then
    echo "  -> No hits below threshold; empty output created: $out_fasta"
    # create empty fasta (or remove existing)
    : > "$out_fasta"
    rm -f "$idfile"
    trap - EXIT
    continue
  fi

  # Extract sequences with headers matching any id (matching first token after '>').
  # Support gzipped protein file: stream via gzip -dc if needed
  if [[ "$PROTEIN" == *.gz ]]; then
    # process stdin (piped)
    eval "$(read_stream_cmd "$PROTEIN")" | \
      awk -v idfile="$idfile" '
        BEGIN {
          while ( (getline line < idfile) > 0 ) {
            ids[line]=1
          }
          close(idfile)
          keep=0
        }
        /^>/ {
          # get first token after '>'
          header=$0
          split(header, a, /[ \t]/)
          id = substr(a[1], 2)
          keep = (id in ids)
        }
        {
          if (keep) print $0
        }
      ' > "$out_fasta"
  else
    awk -v idfile="$idfile" '
      BEGIN {
        while ( (getline line < idfile) > 0 ) {
          ids[line]=1
        }
        close(idfile)
        keep=0
      }
      /^>/ {
        header=$0
        split(header, a, /[ \t]/)
        id = substr(a[1], 2)
        keep = (id in ids)
      }
      {
        if (keep) print $0
      }
    ' "$PROTEIN" > "$out_fasta"
  fi

  # report counts: how many sequences written
  seqcount=$(grep -c '^>' "$out_fasta" || true)
  echo "  -> Extracted $seqcount sequences (found $idcount unique IDs in hits file)."

  rm -f "$idfile"
  trap - EXIT
done

exit 0
