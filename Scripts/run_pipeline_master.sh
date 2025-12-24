#!/usr/bin/env bash
set -euo pipefail

# Ejecución desde el root del repositorio
PIPE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PIPE_ROOT"

# -------- CONFIG --------
CUTINASE_ENV="cutinase_pipe"
SIGNALP_ENV="signalp6"

CUTINASE_PFAM="PF01083"
MAFFT_MODE="auto"
CPU=4

REFS="ref_cutinases.fasta"  # en data/references/
EVALUE_THR="1e-10"

# Summary script (ajusta si tu archivo se llama distinto)
SUMMARY_SCRIPT="scripts/05_summary/01_make_summary.py"

echo "Repo        : $PIPE_ROOT"
echo "cutinase env : $CUTINASE_ENV"
echo "signalp env  : $SIGNALP_ENV"
echo "PFAM id      : $CUTINASE_PFAM"
echo "MAFFT mode   : $MAFFT_MODE"
echo "CPU          : $CPU"
echo "Refs         : data/references/$REFS"
echo "E-value thr  : $EVALUE_THR"
echo

# 1) Verificación - Conda
command -v conda >/dev/null 2>&1 || { echo "ERROR: conda no está en PATH." >&2; exit 1; }

# Scripts esperados
for s in \
  scripts/01_hmmsearch/01_run_hmmsearch.sh \
  scripts/01_hmmsearch/02_filter_extract_fasta.sh \
  scripts/02_signalp/01_run_signalp_cleave.sh \
  scripts/03_pfam/01_runpfam_select_cutinase.sh \
  scripts/04_msa/01_make_msas_per_hit.sh
do
  [ -f "$s" ] || { echo "ERROR: falta el script: $s" >&2; exit 1; }
done

[ -f "$SUMMARY_SCRIPT" ] || { echo "ERROR: falta el summary: $SUMMARY_SCRIPT" >&2; exit 1; }
[ -f "data/references/$REFS" ] || { echo "ERROR: no encuentro refs: data/references/$REFS" >&2; exit 1; }

# 2) Lista de proteomas
shopt -s nullglob
PROTEOMES=(data/proteomes/*.faa data/proteomes/*.fasta data/proteomes/*.fa)
if [ "${#PROTEOMES[@]}" -eq 0 ]; then
  echo "ERROR: No encuentro proteomas en data/proteomes/ (.faa/.fasta/.fa)" >&2
  exit 1
fi

# 3) Modelo HMM
HMM_MODEL=(data/hmm_model/*.hmm)
if [ "${#HMM_MODEL[@]}" -eq 0 ]; then
  echo "ERROR: No encuentro un .hmm en data/hmm_model/." >&2
  exit 1
fi

MODEL="${HMM_MODEL[0]}"
echo "Usando HMM: $MODEL"
echo

# -----------------------
# PIPELINE POR MUESTRA
# -----------------------
for p in "${PROTEOMES[@]}"; do
  base="$(basename "$p")"
  base="${base%.gz}"
  sample="${base%.*}"

  echo "============================================================"
  echo "SAMPLE  : $sample"
  echo "Proteome: $p"
  echo "============================================================"

  # Etapa 01: hmmsearch
  conda run -n "$CUTINASE_ENV" \
    bash scripts/01_hmmsearch/01_run_hmmsearch.sh "$MODEL" "$p"

  # Etapa 01b: filtrar + extraer FASTA
  conda run -n "$CUTINASE_ENV" \
    bash scripts/01_hmmsearch/02_filter_extract_fasta.sh -s "$sample" -t "$EVALUE_THR" -c domain

  # Etapa 02: SignalP + corte (EJECUTAR DIRECTO EN signalp6)
  conda run -n "$SIGNALP_ENV" \
    bash scripts/02_signalp/01_run_signalp_cleave.sh -s "$sample"

  # Etapa 03: Pfam hmmscan + filtrar cutinasa
  conda run -n "$CUTINASE_ENV" \
    bash scripts/03_pfam/01_runpfam_select_cutinase.sh -s "$sample" -p "$CUTINASE_PFAM" -c "$CPU"

  # Etapa 04: MSA por hit
  conda run -n "$CUTINASE_ENV" \
    bash scripts/04_msa/01_make_msas_per_hit.sh -s "$sample" -r "$REFS" -m "$MAFFT_MODE" -c "$CPU"

done

# -----------------------
# RESUMEN FINAL (CSV)
# -----------------------
echo
echo "Generando CSV resumen global..."
conda run -n "$CUTINASE_ENV" \
  python "$SUMMARY_SCRIPT" --pfam "$CUTINASE_PFAM"

echo "LISTO ✅  Revisa: results/summary/cutinase_candidates_summary.csv"
