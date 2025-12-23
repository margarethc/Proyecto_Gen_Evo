#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Script: 02_filter_extract_fasta.sh
# Etapa: 01_hmmsearch
#
# Uso:
#   ./02_filter_extract_fasta.sh -s <sample> [-t 1e-20] [-c domain|full]
#
# Requiere:
#   - domtblout en results/01_hmmsearch/<sample>/<sample>.domtblout
#   - proteoma en data/proteomes/<sample>.faa (o .fasta)
#
# Produce:
#   results/01_hmmsearch/<sample>/
#     ├── <sample>_ids_evalue.txt
#     └── <sample>_hits_filtered.fasta
# ------------------------------------------------------------------------------

set -euo pipefail

THRESHOLD="1e-10"
COLUMN="domain"   # domain=i-Evalue por dominio, full=E-value secuencia completa
SAMPLE=""

usage(){
  cat <<EOF >&2
Uso: $0 -s <sample> [-t THRESHOLD] [-c domain|full]
  -s SAMPLE      Prefijo del proteoma (ej. HP332342k)
  -t THRESHOLD   Umbral e-value (default: $THRESHOLD)
  -c COLUMN      domain o full (default: $COLUMN)
EOF
  exit 2
}

while getopts ":s:t:c:h" opt; do
  case "$opt" in
    s) SAMPLE="$OPTARG" ;;
    t) THRESHOLD="$OPTARG" ;;
    c) COLUMN="$OPTARG" ;;
    h) usage ;;
    \?) echo "Opción inválida: -$OPTARG" >&2; usage ;;
    :)  echo "La opción -$OPTARG requiere un argumento." >&2; usage ;;
  esac
done

[ -z "$SAMPLE" ] && usage

# Root del repo (asumiendo ejecución desde cutinasas_pipeline/)
PIPE_ROOT="$(pwd)"
DOMTBL="${PIPE_ROOT}/results/01_hmmsearch/${SAMPLE}/${SAMPLE}.domtblout"

# Proteoma: intenta .faa primero y luego .fasta/.fa
PROTEOME=""
for ext in faa fasta fa fna; do
  cand="${PIPE_ROOT}/data/proteomes/${SAMPLE}.${ext}"
  if [ -f "$cand" ]; then
    PROTEOME="$cand"
    break
  fi
done

if [ ! -f "$DOMTBL" ]; then
  echo "No se encuentra domtblout: $DOMTBL" >&2
  exit 1
fi
if [ -z "$PROTEOME" ]; then
  echo "No se encuentra proteoma para sample '$SAMPLE' en data/proteomes/ (busqué .faa/.fasta/.fa/.fna)" >&2
  exit 1
fi

if ! command -v seqkit >/dev/null 2>&1; then
  echo "seqkit no encontrado. Instálalo (conda/mamba: seqkit) o dime si prefieres seqtk." >&2
  exit 1
fi

case "$COLUMN" in
  domain) COL_IDX=13 ;;  # domtblout: i-Evalue por dominio (col 13)
  full)   COL_IDX=7  ;;  # domtblout: full sequence E-value (col 7)
  *) echo "COLUMN inválida: $COLUMN (usa 'domain' o 'full')" >&2; exit 2 ;;
esac

OUTDIR="${PIPE_ROOT}/results/01_hmmsearch/${SAMPLE}"
mkdir -p "$OUTDIR"

IDS_FILE="${OUTDIR}/${SAMPLE}_ids_evalue.txt"
OUT_FASTA="${OUTDIR}/${SAMPLE}_hits_filtered.fasta"

echo "==> Filtrando: $DOMTBL"
echo "    sample     : $SAMPLE"
echo "    threshold  : $THRESHOLD"
echo "    columna    : $COLUMN (domtblout col $COL_IDX)"
echo "    proteoma   : $PROTEOME"
echo "    ids        : $IDS_FILE"
echo "    out fasta  : $OUT_FASTA"

# 1) extraer IDs (target_name = col 1 en domtblout)
awk -v col="$COL_IDX" -v thr="$THRESHOLD" '
  BEGIN { thrn = thr + 0 }
  $0 ~ /^#/ { next }
  NF < col { next }
  {
    ev = ($col + 0)
    if (ev < thrn) print $1
  }
' "$DOMTBL" | sort -u > "$IDS_FILE"

idcount=$(wc -l < "$IDS_FILE" | tr -d ' ')
if [ "$idcount" -eq 0 ]; then
  echo "    -> 0 hits bajo el umbral. Creo FASTA vacío: $OUT_FASTA"
  : > "$OUT_FASTA"
  exit 0
fi

# 2) extraer secuencias del proteoma con seqkit (regex anclada)
tmp_re="$(mktemp)"
awk '{print "^"$0"$"}' "$IDS_FILE" > "$tmp_re"
seqkit grep -r -f "$tmp_re" "$PROTEOME" > "$OUT_FASTA"
rm -f "$tmp_re"


seqcount=$(grep -c '^>' "$OUT_FASTA" || true)
echo "    -> IDs únicos: $idcount | Secuencias extraídas: $seqcount"
