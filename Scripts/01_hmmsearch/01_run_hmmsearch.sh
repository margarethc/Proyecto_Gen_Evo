#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Script: 01_run_hmmsearch.sh
# Etapa: 01_hmmsearch
#
# Uso:
#   ./01_run_hmmsearch.sh <modelo.hmm> <proteinas1.faa|fasta> [proteinas2 ...]
#
# Salidas (por proteoma <sample>):
#   results/01_hmmsearch/<sample>/
#     ├── <sample>.domtblout
#     ├── <sample>_hmmsearch.out        
#     └── <sample>_hits.tsv             (tabla legible por dominio)
# ------------------------------------------------------------------------------

set -euo pipefail #si un comando falla o no hay variables el script se detiene

usage() {
  cat <<EOF >&2
Uso: $0 <modelo.hmm> <proteinas.faa|.fasta[.gz]> [proteinas2 ...]
Escribe salidas en: results/01_hmmsearch/<sample>/
EOF
  exit 2
}

if [ "$#" -lt 2 ]; then
  usage
fi

MODEL="$1"
shift

if [ ! -f "$MODEL" ]; then
  echo "Archivo de modelo no encontrado: $MODEL" >&2
  exit 1
fi

if ! command -v hmmsearch >/dev/null 2>&1; then
  echo "hmmsearch no encontrado en PATH. Instala HMMER." >&2
  exit 1
fi

# Root del proyecto = carpeta actual (recomendado ejecutar desde cutinasas_pipeline/)
#Se define la ruta del proyecto y la carpeta de resultados
#Si se ejecuta de otro root, el pipeline se rompe!!!
PIPE_ROOT="$(pwd)"
RESULTS_DIR="${PIPE_ROOT}/results/01_hmmsearch"

mkdir -p "$RESULTS_DIR"

#En caso que uno de los proteomas indicados no exista, el siguiente loop permite su omisión
for PROTEIN in "$@"; do
  if [ ! -f "$PROTEIN" ]; then
    echo "Advertencia: archivo de proteínas no encontrado, se omite: $PROTEIN" >&2
    continue
  fi

# DETERMINAR EL NOMBRE DE LA MUESTRA Y ARCHIVOS (trazabilidad)
  base=$(basename "$PROTEIN") # quita .gz si existiera
  base="${base%.gz}" # quita extensión final (.faa/.fasta/.fa/.fna etc.)
  sample="${base%.*}"

# CREAR UNA CARPETA POR MUESTRA
  outdir="${RESULTS_DIR}/${sample}"
  mkdir -p "$outdir"
  
# NOMBRES DE OUTPUTS
  domtbl="${outdir}/${sample}.domtblout"
  hmmout="${outdir}/${sample}_hmmsearch.out"
  hits_tsv="${outdir}/${sample}_hits.tsv"

  echo "==> hmmsearch: sample='$sample' modelo='$(basename "$MODEL")' proteoma='$base'"
  echo "    domtblout: $domtbl"
  echo "    hits_tsv : $hits_tsv"

  #################
  EJECUTA HMMSEARCH
  #################
  
  if ! hmmsearch --domtblout "$domtbl" "$MODEL" "$PROTEIN" >"$hmmout" 2>&1; then
    echo "hmmsearch falló para $PROTEIN. Revisa: $hmmout" >&2
    continue
  fi

  #################
  Parseo domtblout -> tabla legible por dominio
  #################
  {
    echo -e "query_name\ttarget_name\tdomain_iEvalue\tdomain_bitscore\tfullseq_Evalue\tfullseq_bitscore\thmm_from\thmm_to\tali_from\tali_to\tacc\ttarget_description"
    grep -v '^#' "$domtbl" | awk 'BEGIN{OFS="\t"}{
      desc="";
      if (NF>=23) {
        desc=$23;
        for(i=24;i<=NF;i++) desc=desc " " $i;
      }
      # target_name ($1), query_name ($4), fullseq_Evalue($7), fullseq_bitscore($8),
      # domain_iEvalue($13), domain_bitscore($14), hmm_from($16), hmm_to($17),
      # ali_from($18), ali_to($19), acc($22)
      print $4, $1, $13, $14, $7, $8, $16, $17, $18, $19, $22, desc
    }'
  } > "$hits_tsv"

  echo "    OK: $hits_tsv"
done
