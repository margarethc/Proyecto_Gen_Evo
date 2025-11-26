#!/usr/bin/env bash
# -----------------------------------------------------------------------------------
# Script: filter_and_extract.sh
#
# Uso:
#   ./filter_and_extract.sh [-t THRESHOLD] [-c domain|full] [-o OUTDIR] -- hits1.txt prot1.faa [hits2.txt prot2.faa ...]
#
# Ejemplos:
#   # Umbral por defecto 1e-10, usa e-value por dominio (columna 3)
#   ./filter_and_extract.sh -- sample_hits.txt sample.faa
#
#   # Umbral personalizado y múltiples pares, salida en ./filtered_outputs/
#   ./filter_and_extract.sh -t 1e-34 -c domain -o filtered_outputs -- s1_hits.txt s1.faa s2_hits.txt s2.faa
#
# Notas:
#  - Se espera que los argumentos después de '--' sean pares: archivo de hits + archivo de proteínas.
#  - Los archivos de hits deben estar en formato tabular producido previamente por hmmsearch_tab.sh.
#  - Por defecto se usa la e-value por dominio (columna 3); usar '-c full' para usar e-value por secuencia completa (columna 5).
# -----------------------------------------------------------------------------------

set -euo pipefail  # Terminar en caso de error, uso de variable no definida o fallo en pipeline

# Valores por defecto
THRESHOLD="1e-10"  # Umbral de e-value
COLUMN="domain"     # Columna de e-value a usar: "domain"=col3, "full"=col5
OUTDIR=""
PRINT_HELP=0

# -----------------------------------------------------------------------------------
# Función de ayuda
# -----------------------------------------------------------------------------------
usage() {
  cat <<EOF >&2
Uso: $0 [-t THRESHOLD] [-c domain|full] [-o OUTDIR] -- hits1.txt prot1.faa [hits2.txt prot2.faa ...]
  -t THRESHOLD   Umbral de e-value (default: $THRESHOLD). Se considera hit si evalue < THRESHOLD
  -c COLUMN      Columna a usar: domain (col3) o full (col5). Default: $COLUMN
  -o OUTDIR      Directorio donde se guardan los archivos filtrados (se crea si no existe)
  --             Fin de opciones; luego se indican los pares: hits protein hits protein ...
EOF
  exit 2
}

# -----------------------------------------------------------------------------------
# Parseo de opciones
# -----------------------------------------------------------------------------------
while getopts ":t:c:o:h" opt; do
  case "$opt" in
    t) THRESHOLD="$OPTARG" ;;
    c) COLUMN="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    h) PRINT_HELP=1 ;;
    \?) echo "Opción inválida: -$OPTARG" >&2; usage ;;
    :)  echo "La opción -$OPTARG requiere un argumento." >&2; usage ;;
  esac
done
shift $((OPTIND - 1))

if [ "$PRINT_HELP" -eq 1 ]; then
  usage
fi

# -----------------------------------------------------------------------------------
# Validar que haya al menos un par hits+protein
# -----------------------------------------------------------------------------------
if [ "$#" -lt 2 ]; then
  usage
fi

# Crear directorio de salida si se indicó
if [ -n "$OUTDIR" ]; then
  mkdir -p "$OUTDIR"
fi

# -----------------------------------------------------------------------------------
# Determinar índice de columna para el e-value
# -----------------------------------------------------------------------------------
case "$COLUMN" in
  domain) COL_IDX=3 ;;  # usar columna 3: e-value por dominio
  full)   COL_IDX=5 ;;  # usar columna 5: e-value por secuencia completa
  *)
    echo "Columna inválida: $COLUMN. Elige 'domain' o 'full'." >&2
    exit 2
    ;;
esac

# -----------------------------------------------------------------------------------
# Validar que haya número par de argumentos (pares hits+protein)
# -----------------------------------------------------------------------------------
if [ $(( $# % 2 )) -ne 0 ]; then
  echo "Se esperan pares de argumentos: hits1 prot1 [hits2 prot2 ...]" >&2
  usage
fi

# -----------------------------------------------------------------------------------
# Función helper para leer archivos posiblemente gzipped
# -----------------------------------------------------------------------------------
read_stream_cmd() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    # Usar gzip -dc si está disponible, sino zcat
    if command -v gzip >/dev/null 2>&1; then
      printf "gzip -dc %q" "$f"
    else
      printf "zcat %q" "$f"
    fi
  else
    printf "cat %q" "$f"
  fi
}

# -----------------------------------------------------------------------------------
# Loop sobre pares de archivos: hits + protein
# -----------------------------------------------------------------------------------
while [ "$#" -gt 0 ]; do
  HITFILE="$1"; shift
  PROTEIN="$1"; shift

  # Validación de existencia de archivos
  if [ ! -f "$HITFILE" ]; then
    echo "Advertencia: archivo de hits no encontrado, se omite: $HITFILE" >&2
    continue
  fi
  if [ ! -f "$PROTEIN" ]; then
    echo "Advertencia: archivo de proteínas no encontrado, se omite: $PROTEIN" >&2
    continue
  fi

  # ---------------------------------------------------------------------------------
  # Determinar nombre de salida FASTA filtrado
  # ---------------------------------------------------------------------------------
  base=$(basename "$PROTEIN")
  [[ "$base" == *.gz ]] && base="${base%.gz}"   # eliminar extensión .gz
  name="${base%.*}"                             # eliminar extensión final (.faa/.fasta)
  out_fasta="${name}_filtered.fasta"
  [[ -n "$OUTDIR" ]] && out_fasta="${OUTDIR%/}/${out_fasta}"

  echo "Procesando par:"
  echo "  hits:    $HITFILE"
  echo "  protein: $PROTEIN"
  echo "  threshold: $THRESHOLD (usando columna $COL_IDX - $COLUMN)"
  echo "  salida:  $out_fasta"

  # ---------------------------------------------------------------------------------
  # Extraer IDs de secuencias que cumplen el umbral de e-value
  # ---------------------------------------------------------------------------------
  idfile="$(mktemp)"  # archivo temporal para almacenar IDs
  trap 'rm -f "$idfile"' EXIT

  eval "$(read_stream_cmd "$HITFILE")" | \
    awk -v col="$COL_IDX" -v thr="$THRESHOLD" '
      BEGIN{
        thrn = thr + 0  # convertir e-value a numérico (soporta notación científica)
      }
      /^#/ { next }                             # saltar comentarios HMMER
      NR==1 && tolower($0) ~ /query_name/ { next } # saltar encabezado si existe
      NF < col { next }                         # saltar líneas mal formateadas
      {
        ev = ($col + 0)
        if (ev < thrn) {
          print $2   # columna 2 = target_name
        }
      }
    ' | sort -u > "$idfile"

  # Contar cuántos IDs válidos se encontraron
  idcount=$(wc -l < "$idfile" | tr -d ' ')
  if [ "$idcount" -eq 0 ]; then
    echo "  -> No se encontraron hits por debajo del umbral; salida vacía: $out_fasta"
    : > "$out_fasta"  # crear archivo vacío
    rm -f "$idfile"
    trap - EXIT
    continue
  fi

  # ---------------------------------------------------------------------------------
  # Extraer secuencias correspondientes de archivo FASTA (soporta gz)
  # ---------------------------------------------------------------------------------
  if [[ "$PROTEIN" == *.gz ]]; then
    eval "$(read_stream_cmd "$PROTEIN")" | \
      awk -v idfile="$idfile" '
        BEGIN {
          while ((getline line < idfile) > 0) ids[line]=1
          close(idfile)
          keep=0
        }
        /^>/ {
          header=$0
          split(header, a, /[ \t]/)
          id = substr(a[1], 2)  # primer token del header
          keep = (id in ids)
        }
        { if (keep) print $0 }
      ' > "$out_fasta"
  else
    awk -v idfile="$idfile" '
      BEGIN {
        while ((getline line < idfile) > 0) ids[line]=1
        close(idfile)
        keep=0
      }
      /^>/ {
        header=$0
        split(header, a, /[ \t]/)
        id = substr(a[1], 2)
        keep = (id in ids)
      }
      { if (keep) print $0 }
    ' "$PROTEIN" > "$out_fasta"
  fi

  # ---------------------------------------------------------------------------------
  # Reportar número de secuencias extraídas
  # ---------------------------------------------------------------------------------
  seqcount=$(grep -c '^>' "$out_fasta" || true)
  echo "  -> Se extrajeron $seqcount secuencias (encontradas $idcount IDs únicos en hits)."

  rm -f "$idfile"
  trap - EXIT
done

exit 0
