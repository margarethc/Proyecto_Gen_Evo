#!/usr/bin/env bash

####################################################################################################
# Run the script from inside GE_cut/:
#
#       a) Process **all** FASTA files in filtered_outputs/:
#
#              ./run_signalp_cleave.sh
#
#       b) Process **only one** specific file:
#
#              ./run_signalp_cleave.sh CAL9_filtered.fasta
#
#       c) Process **multiple selected** files:
#
#              ./run_signalp_cleave.sh CAL9_filtered.fasta RUE08_filtered.faa FUSOX_filtered.fa
#
#

####################################################################################################

set -euo pipefail

# Directory with your filtered FASTA files
INPUT_DIR="filtered_outputs"

# Base directory for SignalP outputs
RESULTS_BASE="signalp_results"
mkdir -p "$RESULTS_BASE"

# Determine which files to process:
# - No arguments  → process all FASTA files in filtered_outputs/
# - One or more   → process only the specified files (by name, inside filtered_outputs/)
if [ "$#" -eq 0 ]; then
    # No arguments: process all FASTA-like files
    FILES_TO_PROCESS=("$INPUT_DIR"/*.fasta "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.faa)
else
    # User provided one or more filenames
    FILES_TO_PROCESS=()
    for arg in "$@"; do
        fa="$INPUT_DIR/$arg"
        if [ ! -f "$fa" ]; then
            echo "ERROR: File '$arg' not found in $INPUT_DIR/"
            exit 1
        fi
        FILES_TO_PROCESS+=("$fa")
    done
fi

processed_any=false

for fa in "${FILES_TO_PROCESS[@]}"; do
    # Handle the case when the glob didn't match anything in "all files" mode
    if [ ! -e "$fa" ]; then
        continue
    fi

    processed_any=true

    fname=$(basename "$fa")            # e.g. CAL9_filtered.fasta
    base_noext="${fname%%.*}"         # e.g. CAL9_filtered
    core="${base_noext%_filtered}"    # e.g. CAL9  (removes "_filtered" if present)

    echo "==============================================="
    echo "Processing file: $fa"
    echo "  Sample name:   $core"

    # Per-sample output directory for SignalP raw outputs
    outdir="$RESULTS_BASE/$core"
    mkdir -p "$outdir"

    ########################################
    # 1) Run SignalP 6 (slow-sequential)
    ########################################
    signalp6 \
        --fastafile "$fa" \
        --organism euk \
        --output_dir "$outdir" \
        --format none \
        --mode slow-sequential

    gff="$outdir/output.gff3"

    if [ ! -f "$gff" ]; then
        echo "  WARNING: No output.gff3 found for $fa, skipping cleavage step."
        continue
    fi

    ########################################
    # 2) Extract SP coordinates from GFF3
    ########################################
    tmp_table=$(mktemp)

    # GFF3: seqid (col1), feature type (col3), start (col4), end (col5)
    # We keep only feature type "signal_peptide"
    awk -F'\t' '
        $0 !~ /^#/ && $3 == "signal_peptide" {
            print $1 "\t" $4 "\t" $5
        }
    ' "$gff" > "$tmp_table"

    # If no signal peptides found, skip this file
    if [ ! -s "$tmp_table" ]; then
        echo "  No signal peptides detected in $fa."
        rm -f "$tmp_table"
        continue
    fi

    ########################################
    # 3) Create cleaved FASTA + per-sample summary
    ########################################
    # Output FASTA: e.g. CAL9_filtered_signalp.fasta
    cleaved_fasta="${base_noext}_signalp.fasta"

    # Summary file: e.g. CAL9_signalP_summary.tsv
    summary_file="${core}_signalP_summary.tsv"

    python <<PY
import textwrap

fasta_path = "$fa"
table_path = "$tmp_table"
out_fasta_path = "$cleaved_fasta"
summary_path = "$summary_file"

# Load cleavage positions
# table has: seqid \t sp_start \t sp_end
cleavage = {}
with open(table_path) as t:
    for line in t:
        line = line.strip()
        if not line:
            continue
        seqid, sp_start, sp_end = line.split("\t")
        # If there are multiple SP entries for the same seqid, keep the first
        if seqid not in cleavage:
            cleavage[seqid] = (int(sp_start), int(sp_end))

def read_fasta(path):
    header_id = None
    header_full = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header_id is not None:
                    yield header_id, header_full, "".join(seq_lines)
                header_full = line[1:]              # full header (without '>')
                header_id = header_full.split()[0]  # ID = first token
                seq_lines = []
            else:
                seq_lines.append(line)
    if header_id is not None:
        yield header_id, header_full, "".join(seq_lines)

with open(out_fasta_path, "w") as out_fa, open(summary_path, "w") as out_sum:
    out_sum.write("seqid\tsp_start\tsp_end\tcleavage_after_aa\toriginal_len\tnew_len\n")

    for sid, full_header, seq in read_fasta(fasta_path):
        if sid not in cleavage:
            # No signal peptide predicted for this sequence
            continue

        sp_start, sp_end = cleavage[sid]
        # Cleavage is AFTER the last AA of the signal peptide
        cut_pos = sp_end  # 1-based index
        # Python string is 0-based, so new sequence starts at cut_pos
        new_seq = seq[cut_pos:]

        if len(new_seq) == 0:
            # Resulting sequence empty → skip
            continue

        # Write cleaved sequence to FASTA, keep original header but mark as cleaved
        out_fa.write(f">{full_header} | signalp_cleaved\n")
        for chunk in textwrap.wrap(new_seq, 60):
            out_fa.write(chunk + "\n")

        # Write summary line
        out_sum.write(
            f"{sid}\t{sp_start}\t{sp_end}\t{cut_pos}\t{len(seq)}\t{len(new_seq)}\n"
        )
PY

    rm -f "$tmp_table"

    echo "  -> Cleaved FASTA: $cleaved_fasta"
    echo "  -> Summary file:  $summary_file"
    echo "  -> Raw SignalP outputs in: $outdir"
done

if [ "$processed_any" = false ]; then
    echo "No FASTA files were processed. Check that filtered_outputs/ contains .fasta/.fa/.faa files,"
    echo "or that the filenames you provided as arguments exist in that folder."
fi

echo "Done."
