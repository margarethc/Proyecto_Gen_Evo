#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Script: 01_run_pfam_select_cutinase.sh
# Etapa: 03_pfam
#
# Uso:
#   ./01_run_pfam_select_cutinase.sh -s <sample> -p <PFAM_ID> [-c 4]
#
# Entrada:
#   results/02_signalp/<sample>/<sample>_signalp_trimmed.fasta
#
# Salidas:
#   results/03_pfam/<sample>/
#     ├── <sample>_pfam.domtblout
#     ├── <sample>_pfam.hmmscan.out
#     ├── <sample>_pfam_cutinase_ids.txt
#     └── <sample>_pfam_filtered.fasta
# ------------------------------------------------------------------------------

set -euo pipefail

SAMPLE=""
PFAM_ID=""
CPU=4

usage(){
  cat <<EOF >&2
Uso: $0 -s <sample> -p <PFAM_ID> [-c CPU]
  -s SAMPLE    Prefijo del proteoma (ej. HP332342k)
  -p PFAM_ID   Accession Pfam del dominio cutinasa (ej. PF01083)
  -c CPU       CPUs para hmmscan (default: $CPU)
EOF
  exit 2
}

while getopts ":s:p:c:h" opt; do
  case "$opt" in
    s) SAMPLE="$OPTARG" ;;
    p) PFAM_ID="$OPTARG" ;;
    c) CPU="$OPTARG" ;;
    h) usage ;;
    \?) echo "Opción inválida: -$OPTARG" >&2; usage ;;
    :)  echo "La opción -$OPTARG requiere un argumento." >&2; usage ;;
  esac
done

[ -z "$SAMPLE" ] && usage
[ -z "$PFAM_ID" ] && usage

PIPE_ROOT="$(pwd)"

# Entrada: SOLO secretadas (con péptido señal)
IN_FASTA="${PIPE_ROOT}/results/02_signalp/${SAMPLE}/${SAMPLE}_signalp_trimmed.fasta"
if [ ! -f "$IN_FASTA" ]; then
  echo "ERROR: No encuentro el FASTA secretado: $IN_FASTA" >&2
  exit 1
fi

# Uso de la base de datos de PFam
PFAM_HMM="${PIPE_ROOT}/databases/pfam/Pfam-A.hmm"
if [ ! -f "$PFAM_HMM" ]; then
  # fallback por si lo tienes directo en databases/
  PFAM_HMM="${PIPE_ROOT}/databases/Pfam-A.hmm"
fi
if [ ! -f "$PFAM_HMM" ]; then
  echo "ERROR: No encuentro Pfam-A.hmm en databases/pfam/ ni databases/." >&2
  exit 1
fi

if ! command -v hmmscan >/dev/null 2>&1; then
  echo "ERROR: hmmscan no está en PATH!!" >&2
  exit 1
fi

if ! command -v seqkit >/dev/null 2>&1; then
  echo "ERROR: seqkit no está en PATH!!" >&2
  exit 1
fi

#Se define el directorio para los outputs
OUTDIR="${PIPE_ROOT}/results/03_pfam/${SAMPLE}"
mkdir -p "$OUTDIR"

#Se define los nombres de los outputs
DOMTBL="${OUTDIR}/${SAMPLE}_pfam.domtblout"
HMMOUT="${OUTDIR}/${SAMPLE}_pfam.hmmscan.out"
IDS="${OUTDIR}/${SAMPLE}_pfam_cutinase_ids.txt"
OUT_FASTA="${OUTDIR}/${SAMPLE}_pfam_filtered.fasta"

echo "==> Pfam hmmscan + selección dominio"
echo "    sample   : $SAMPLE"
echo "    input    : $IN_FASTA"
echo "    pfam_db  : $PFAM_HMM"
echo "    pfam_id  : $PFAM_ID"
echo "    outdir   : $OUTDIR"

# 1) hmmscan contra Pfam-A
hmmscan --cpu "$CPU" \
  --domtblout "$DOMTBL" \
  "$PFAM_HMM" "$IN_FASTA" > "$HMMOUT"

# 2) filtrar hits por Pfam accession
# domtblout: $1=target name, $2=target accession, $4=query name
awk -v pf="$PFAM_ID" '
  $0 ~ /^#/ { next }
  {
    acc=$2
    sub(/\..*/, "", acc)   # PF01083.23 -> PF01083 # Dominio para cutinasas!!
    if (acc == pf) print $4
  }
' "$DOMTBL" | sort -u > "$IDS"

idcount=$(wc -l < "$IDS" | tr -d ' ')
if [ "$idcount" -eq 0 ]; then
  echo "    -> 0 secuencias con dominio $PFAM_ID. Creo FASTA vacío: $OUT_FASTA"
  : > "$OUT_FASTA"
  exit 0
fi

# 3) se  extrae del FASTA secretado las secuencias que tienen el dominio
tmp_re="$(mktemp)"
awk '{print "^"$0"$"}' "$IDS" > "$tmp_re"
seqkit grep -r -f "$tmp_re" "$IN_FASTA" > "$OUT_FASTA"
rm -f "$tmp_re"

# Para el reporte final
seqcount=$(grep -c '^>' "$OUT_FASTA" || true)
echo "    -> IDs con $PFAM_ID: $idcount | Secuencias en FASTA final: $seqcount"
echo "    OK: $OUT_FASTA"
