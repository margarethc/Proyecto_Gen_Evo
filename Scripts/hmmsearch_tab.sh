#!/usr/bin/env bash
# -----------------------------------------------------------------------------------
# Script: hmmsearch_tab.sh
# Autor: (Tu nombre)
#
# Descripción:
#   Este script ejecuta HMMER (hmmsearch) usando un perfil HMM (*.hmm) contra uno o más
#   archivos de proteínas (*.faa / *.fasta) y genera un archivo tabulado (*_hits.txt)
#   por cada archivo de proteínas. Cada fila corresponde a un dominio detectado.
#
# Uso:
#   ./hmmsearch_tab.sh modelo.hmm prot1.fasta prot2.faa ...
#   ./hmmsearch_tab.sh modelo.hmm prot1.fasta prot2.faa /ruta/a/salida/   # si el último argumento es un directorio existente, se usa como carpeta de salida
#
# Notas:
#  - Soporta archivos .gz (elimina la extensión .gz al nombrar la salida)
#  - Requiere que hmmsearch esté en el PATH
#  - Produce salida por dominio usando --domtblout (formato TSV)
# -----------------------------------------------------------------------------------

set -euo pipefail  # Termina si hay errores, si se usan variables no definidas, o falla un comando en pipeline

# -----------------------------------------------------------------------------------
# Función de uso
# -----------------------------------------------------------------------------------
usage() {
  cat <<EOF >&2
Uso: $0 <modelo.hmm> <proteinas.faa|.fasta[.gz]> [proteinas2 ...] [out_dir]
Si el último argumento es un directorio existente se usará como carpeta de salida.
EOF
  exit 2
}

# Comprobación mínima de argumentos
if [ "$#" -lt 2 ]; then
  usage
fi

# Primer argumento = modelo HMM
MODEL="$1"
shift

# Verifica que el archivo HMM exista
if [ ! -f "$MODEL" ]; then
  echo "Archivo de modelo no encontrado: $MODEL" >&2
  exit 1
fi

# -----------------------------------------------------------------------------------
# Determinar si el último argumento es un directorio de salida
# -----------------------------------------------------------------------------------
last_arg="${!#}"  # último argumento
OUTDIR=""
if [ -d "$last_arg" ]; then
  OUTDIR="${last_arg%/}"  # elimina slash final si existe
  # Elimina el último argumento de la lista de archivos
  set -- "${@:1:$(($#-1))}"
fi

# Lista de archivos de proteínas restantes
if [ "$#" -lt 1 ]; then
  usage
fi

# -----------------------------------------------------------------------------------
# Verifica que hmmsearch esté disponible
# -----------------------------------------------------------------------------------
if ! command -v hmmsearch >/dev/null 2>&1; then
  echo "hmmsearch no encontrado en PATH. Instala HMMER (http://hmmer.org/)" >&2
  exit 1
fi

# -----------------------------------------------------------------------------------
# Función auxiliar para calcular nombre de salida por cada archivo de proteínas
# -----------------------------------------------------------------------------------
default_out_from_protein() {
  local p="$1"
  local base tmpname name

  base=$(basename "$p")  # Extrae nombre de archivo sin ruta
  tmpname="$base"

  # Si está comprimido (.gz), eliminar extensión
  if [[ "$tmpname" == *.gz ]]; then
    tmpname="${tmpname%.gz}"
  fi

  # Quita la extensión final (.fasta, .faa)
  name="${tmpname%.*}"

  # Construye ruta completa de salida si hay OUTDIR
  if [ -n "$OUTDIR" ]; then
    echo "${OUTDIR}/${name}_hits.txt"
  else
    echo "${name}_hits.txt"
  fi
}

# -----------------------------------------------------------------------------------
# Procesa cada archivo de proteínas
# -----------------------------------------------------------------------------------
for PROTEIN in "$@"; do
  # Salta archivos que no existan
  if [ ! -f "$PROTEIN" ]; then
    echo "Advertencia: archivo de proteínas no encontrado, se omite: $PROTEIN" >&2
    continue
  fi

  # Nombre de salida y archivo temporal para domtblout
  OUT="$(default_out_from_protein "$PROTEIN")"
  TMP_DOM=$(mktemp --suffix=.domtblout)

  # Función de limpieza automática del archivo temporal
  cleanup_file() {
    rm -f "$TMP_DOM"
  }
  trap cleanup_file EXIT

  echo "Ejecutando hmmsearch: modelo='$MODEL' proteínas='$PROTEIN' -> salida='$OUT'"

  # Ejecuta hmmsearch, guarda domtblout en TMP_DOM, descarta stdout/stderr
  if ! hmmsearch --domtblout "$TMP_DOM" "$MODEL" "$PROTEIN" >/dev/null 2>/dev/null; then
    echo "hmmsearch falló para $PROTEIN. Se omite." >&2
    rm -f "$TMP_DOM"
    trap - EXIT
    continue
  fi

  # ---------------------------------------------------------------------------------
  # Parseo del domtblout para generar tabla legible por dominio
  # Columnas de salida: query_name, target_name, domain_iEvalue, domain_bitscore,
  #                     fullseq_Evalue, fullseq_bitscore, hmm_from, hmm_to,
  #                     ali_from, ali_to, acc, target_description
  # ---------------------------------------------------------------------------------
  {
    echo -e "query_name\ttarget_name\tdomain_iEvalue\tdomain_bitscore\tfullseq_Evalue\tfullseq_bitscore\thmm_from\thmm_to\tali_from\tali_to\tacc\ttarget_description"
    grep -v '^#' "$TMP_DOM" | awk 'BEGIN{OFS="\t"} {
      desc="";
      if (NF>=23) {
        desc=$23;
        for(i=24;i<=NF;i++) desc=desc " " $i;
      }
      # Columnas HMMER domtblout: target_name ($1), target_accession ($2), tlen ($3), query_name ($4), etc.
      print $4, $1, $13, $14, $7, $8, $16, $17, $18, $19, $22, desc
    }'
  } > "$OUT"

  echo "Resultados tabulados escritos en: $OUT"

  # Limpieza de archivo temporal
  rm -f "$TMP_DOM"
  trap - EXIT
done
