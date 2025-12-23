#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Script: 01_run_signalp_cleave.sh
# Etapa: 02_signalp
#
# Uso:
#   ./01_run_signalp_cleave.sh -s HP332342k
#   ./01_run_signalp_cleave.sh -s HP332342k -k   # (opcional) mantiene no-SP sin recortar
#
# Entrada esperada:
#   results/01_hmmsearch/<sample>/<sample>_hits_filtered.fasta
#
# Salidas:
#   results/02_signalp/<sample>/
#     ├── <sample>_signalp_out/                 (salida cruda de SignalP)
#     ├── <sample>_signalP_summary.tsv          (tabla SP_start/end)
#     ├── <sample>_signalp_trimmed.fasta        (solo secuencias con SP, recortadas)
#     └── <sample>_signalp_noSP.fasta           (solo secuencias sin SP)  [opcional]
#     └── <sample>_signalp_kept.fasta           (SP recortadas + noSP sin cambio) [si -k]
# ------------------------------------------------------------------------------

set -euo pipefail

SAMPLE=""
KEEP_NOSP=0

usage(){
  cat <<EOF >&2
Uso: $0 -s <sample> [-k]
  -s SAMPLE   Prefijo del proteoma (ej. HP332342k)
  -k          Mantener también las secuencias sin péptido señal (sin recortar) en *_signalp_kept.fasta
EOF
  exit 2
}

while getopts ":s:kh" opt; do
  case "$opt" in
    s) SAMPLE="$OPTARG" ;;
    k) KEEP_NOSP=1 ;;
    h) usage ;;
    \?) echo "Opción inválida: -$OPTARG" >&2; usage ;;
    :)  echo "La opción -$OPTARG requiere un argumento." >&2; usage ;;
  esac
done

[ -z "$SAMPLE" ] && usage

PIPE_ROOT="$(pwd)"

IN_FASTA="${PIPE_ROOT}/results/01_hmmsearch/${SAMPLE}/${SAMPLE}_hits_filtered.fasta"
if [ ! -f "$IN_FASTA" ]; then
  echo "ERROR: No encuentro el FASTA filtrado: $IN_FASTA" >&2
  echo "       (Debe existir desde la etapa 01)" >&2
  exit 1
fi

OUTDIR="${PIPE_ROOT}/results/02_signalp/${SAMPLE}"
RAW_DIR="${OUTDIR}/${SAMPLE}_signalp_out"
mkdir -p "$RAW_DIR"

SUMMARY_TSV="${OUTDIR}/${SAMPLE}_signalP_summary.tsv"
TRIMMED_FASTA="${OUTDIR}/${SAMPLE}_signalp_trimmed.fasta"
NOSP_FASTA="${OUTDIR}/${SAMPLE}_signalp_noSP.fasta"
KEPT_FASTA="${OUTDIR}/${SAMPLE}_signalp_kept.fasta"

echo "==> SignalP6 + corte SP"
echo "    sample      : $SAMPLE"
echo "    input fasta : $IN_FASTA"
echo "    outdir      : $OUTDIR"
echo "    raw dir     : $RAW_DIR"

# Ejecutar SignalP6 desde el env dedicado (recomendado)
# Ajustes “seguros” para tu WSL (RAM ~7.4GiB):
SIGNALP_MODELS="/home/margaret/Proy_gen_evo/signalp6_slow_sequential/signalp-6-package/models"

conda run -n signalp6 signalp6 \
  --fastafile "$IN_FASTA" \
  --output_dir "$RAW_DIR" \
  --organism eukarya \
  --mode slow-sequential \
  --format none \
  --model_dir "$SIGNALP_MODELS" \
  --bsize 2 \
  --write_procs 2 \
  --torch_num_threads 4

GFF="${RAW_DIR}/output.gff3"
if [ ! -f "$GFF" ]; then
  echo "ERROR: No se generó $GFF. Revisa el contenido de: $RAW_DIR" >&2
  exit 1
fi

TMP_TABLE="$(mktemp)"
trap 'rm -f "$TMP_TABLE"' EXIT

# Extraer coordenadas de signal_peptide (seqid, start, end)
awk -F'\t' '
  $0 !~ /^#/ && $3 == "signal_peptide" {
    split($1,a," ");
    print a[1] "\t" $4 "\t" $5
  }
' "$GFF" > "$TMP_TABLE"

# Si no hay SP, igual generamos outputs (vacíos) para trazabilidad
if [ ! -s "$TMP_TABLE" ]; then
  echo "    -> No se detectaron péptidos señal en $SAMPLE."
  : > "$TRIMMED_FASTA"
  : > "$SUMMARY_TSV"
  if [ "$KEEP_NOSP" -eq 1 ]; then
    cp "$IN_FASTA" "$KEPT_FASTA"
  fi
  exit 0
fi

# Cortar FASTA y generar resumen
python - <<PY
import textwrap

fasta_path = "$IN_FASTA"
table_path = "$TMP_TABLE"
trimmed_path = "$TRIMMED_FASTA"
nosp_path = "$NOSP_FASTA"
kept_path = "$KEPT_FASTA"
summary_path = "$SUMMARY_TSV"
keep_nosp = bool($KEEP_NOSP)

# cargar posiciones SP (primera ocurrencia por seqid)
cleavage = {}
with open(table_path) as t:
    for line in t:
        line = line.strip()
        if not line:
            continue
        seqid, sp_start, sp_end = line.split("\\t")
        if seqid not in cleavage:
            cleavage[seqid] = (int(sp_start), int(sp_end))

def read_fasta(path):
    header_id = None
    header_full = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\\n")
            if not line:
                continue
            if line.startswith(">"):
                if header_id is not None:
                    yield header_id, header_full, "".join(seq_lines)
                header_full = line[1:]
                header_id = header_full.split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header_id is not None:
            yield header_id, header_full, "".join(seq_lines)

def write_fasta_record(handle, header, seq):
    handle.write(f">{header}\\n")
    for chunk in textwrap.wrap(seq, 60):
        handle.write(chunk + "\\n")

n_trimmed = 0
n_nosp = 0

with open(trimmed_path, "w") as out_trim, \
     open(nosp_path, "w") as out_nosp, \
     open(summary_path, "w") as out_sum:

    out_sum.write("seqid\\tsp_start\\tsp_end\\tcleavage_after_aa\\toriginal_len\\tnew_len\\n")

    kept_handle = open(kept_path, "w") if keep_nosp else None

    for sid, full_header, seq in read_fasta(fasta_path):
        if sid in cleavage:
            sp_start, sp_end = cleavage[sid]
            cut_pos = sp_end  # 1-based
            new_seq = seq[cut_pos:]  # python slice usa 0-based; cut_pos=sp_end corta después del SP

            if len(new_seq) == 0:
                continue

            write_fasta_record(out_trim, full_header + " | signalp_cleaved", new_seq)
            out_sum.write(f"{sid}\\t{sp_start}\\t{sp_end}\\t{cut_pos}\\t{len(seq)}\\t{len(new_seq)}\\n")
            n_trimmed += 1

            if kept_handle is not None:
                write_fasta_record(kept_handle, full_header + " | signalp_cleaved", new_seq)
        else:
            n_nosp += 1
            write_fasta_record(out_nosp, full_header + " | no_signal_peptide", seq)
            if kept_handle is not None:
                write_fasta_record(kept_handle, full_header + " | no_signal_peptide", seq)

    if kept_handle is not None:
        kept_handle.close()

print(f"[SignalP] trimmed (SP): {n_trimmed} | noSP: {n_nosp}")
PY

echo "    OK:"
echo "      - SP recortadas : $TRIMMED_FASTA"
echo "      - noSP         : $NOSP_FASTA"
echo "      - resumen      : $SUMMARY_TSV"
echo "      - raw          : $RAW_DIR"
if [ "$KEEP_NOSP" -eq 1 ]; then
  echo "      - kept (SP+noSP): $KEPT_FASTA"
fi
