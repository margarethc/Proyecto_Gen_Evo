#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Script: 01_make_msas_per_hit.sh
# Etapa: 04_msa
#
# Uso:
#   ./01_make_msas_per_hit.sh -s <sample> -r <refs.fasta> [-m auto|hq] [-c CPU]
#
# Entrada:
#   results/03_pfam/<sample>/<sample>_pfam_filtered.fasta
#   data/references/<refs.fasta>
#
# Salida:
#   results/04_msa/<sample>/
#     ├── <sample>_query_plus_refs/<hit>.fasta
#     └── <sample>_alignments/<hit>_aln.fas
# ------------------------------------------------------------------------------

set -euo pipefail

SAMPLE=""
REFS=""
MODE="auto"   # auto o hq
CPU=4

usage(){
  cat <<EOF >&2
Uso: $0 -s <sample> -r <refs.fasta> [-m auto|hq] [-c CPU]
  -s SAMPLE     Prefijo del proteoma (ej. HP332342k)
  -r REFS       FASTA de referencias (ruta o nombre dentro de data/references/)
  -m MODE       auto (mafft --auto) o hq (mafft --localpair --maxiterate 1000). Default: $MODE
  -c CPU        CPUs para MAFFT (default: $CPU)
EOF
  exit 2
}

while getopts ":s:r:m:c:h" opt; do
  case "$opt" in
    s) SAMPLE="$OPTARG" ;;
    r) REFS="$OPTARG" ;;
    m) MODE="$OPTARG" ;;
    c) CPU="$OPTARG" ;;
    h) usage ;;
    \?) echo "Opción inválida: -$OPTARG" >&2; usage ;;
    :)  echo "La opción -$OPTARG requiere un argumento." >&2; usage ;;
  esac
done

[ -z "$SAMPLE" ] && usage
[ -z "$REFS" ] && usage

PIPE_ROOT="$(pwd)"

IN_FASTA="${PIPE_ROOT}/results/03_pfam/${SAMPLE}/${SAMPLE}_pfam_filtered.fasta"
if [ ! -f "$IN_FASTA" ]; then
  echo "ERROR: No encuentro el FASTA Pfam filtrado: $IN_FASTA" >&2
  exit 1
fi

# Resolver refs: si no es ruta absoluta, asumir data/references/
if [ -f "$REFS" ]; then
  REFS_FASTA="$REFS"
else
  REFS_FASTA="${PIPE_ROOT}/data/references/${REFS}"
fi

if [ ! -f "$REFS_FASTA" ]; then
  echo "ERROR: No encuentro referencias: $REFS_FASTA" >&2
  exit 1
fi

if ! command -v mafft >/dev/null 2>&1; then
  echo "ERROR: mafft no está en PATH." >&2
  exit 1
fi
if ! command -v seqkit >/dev/null 2>&1; then
  echo "ERROR: seqkit no está en PATH. Instala seqkit en cutinase_pipe." >&2
  exit 1
fi

OUTDIR="${PIPE_ROOT}/results/04_msa/${SAMPLE}"
QDIR="${OUTDIR}/${SAMPLE}_query_plus_refs"
ADIR="${OUTDIR}/${SAMPLE}_alignments"
mkdir -p "$QDIR" "$ADIR"

echo "==> MSA por hit + referencias"
echo "    sample   : $SAMPLE"
echo "    input    : $IN_FASTA"
echo "    refs     : $REFS_FASTA"
echo "    mode     : $MODE"
echo "    outdir   : $OUTDIR"

# Obtener IDs del FASTA final (hits)
mapfile -t IDS < <(seqkit seq -n "$IN_FASTA")

if [ "${#IDS[@]}" -eq 0 ]; then
  echo "    -> No hay secuencias en $IN_FASTA. Nada que alinear."
  exit 0
fi

# Elegir parámetros MAFFT
if [ "$MODE" = "hq" ]; then
  MAFFT_ARGS=(--localpair --maxiterate 1000 --thread "$CPU")
else
  MAFFT_ARGS=(--auto --thread "$CPU")
fi

for id in "${IDS[@]}"; do
  # ID seguro para nombres de archivo
  safe_id="$(echo "$id" | sed 's/[^A-Za-z0-9._-]/_/g')"

  qfa="${QDIR}/${safe_id}.fasta"
  aln="${ADIR}/${safe_id}_aln.fas"

  # Extraer la secuencia de la query de forma robusta
  tmp_id="${QDIR}/.tmp_id_${safe_id}.txt"
  printf "%s\n" "$id" > "$tmp_id"
  seqkit grep -n -f "$tmp_id" "$IN_FASTA" > "$qfa"
  rm -f "$tmp_id"

  # Añadir referencias
  cat "$REFS_FASTA" >> "$qfa"

  # Alinear
  mafft "${MAFFT_ARGS[@]}" "$qfa" > "$aln"

  echo "    OK: $aln"
done
