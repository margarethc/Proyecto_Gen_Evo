#!/usr/bin/env bash

####################################################################################################
# Script para ejecutar SignalP 6 sobre archivos FASTA de proteínas y extraer secuencias
# después del péptido señal predicho (cleavage).  
#
# Uso desde el directorio GE_cut/:
#   a) Procesar todos los FASTA en filtered_outputs/:
#        ./run_signalp_cleave.sh
#   b) Procesar un archivo específico:
#        ./run_signalp_cleave.sh CAL9_filtered.fasta
#   c) Procesar múltiples archivos seleccionados:
#        ./run_signalp_cleave.sh CAL9_filtered.fasta RUE08_filtered.faa FUSOX_filtered.fa
####################################################################################################

set -euo pipefail  # Terminar si hay error, variables no definidas o errores en pipes

# Carpeta que contiene los FASTA filtrados
INPUT_DIR="filtered_outputs"

# Carpeta base donde se guardarán los resultados de SignalP
RESULTS_BASE="signalp_results"
mkdir -p "$RESULTS_BASE"

# -----------------------------
# Determinar qué archivos procesar
# -----------------------------
# Sin argumentos -> todos los FASTA en filtered_outputs/
# Con argumentos -> solo los archivos especificados
if [ "$#" -eq 0 ]; then
    FILES_TO_PROCESS=("$INPUT_DIR"/*.fasta "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.faa)
else
    FILES_TO_PROCESS=()
    for arg in "$@"; do
        fa="$INPUT_DIR/$arg"
        if [ ! -f "$fa" ]; then
            echo "ERROR: Archivo '$arg' no encontrado en $INPUT_DIR/"
            exit 1
        fi
        FILES_TO_PROCESS+=("$fa")
    done
fi

processed_any=false  # Bandera para saber si se procesó al menos un archivo

# -----------------------------
# Bucle principal: procesar cada archivo FASTA
# -----------------------------
for fa in "${FILES_TO_PROCESS[@]}"; do
    # Saltar si el glob no encontró coincidencias
    if [ ! -e "$fa" ]; then
        continue
    fi

    processed_any=true

    fname=$(basename "$fa")            # Nombre del archivo, ej: CAL9_filtered.fasta
    base_noext="${fname%%.*}"         # Sin extensión, ej: CAL9_filtered
    core="${base_noext%_filtered}"    # Quita "_filtered", ej: CAL9

    echo "==============================================="
    echo "Procesando archivo: $fa"
    echo "  Nombre de muestra: $core"

    # Carpeta de salida específica para cada muestra
    outdir="$RESULTS_BASE/$core"
    mkdir -p "$outdir"

    ########################################
    # 1) Ejecutar SignalP 6 (modo slow-sequential)
    ########################################
    signalp6 \
        --fastafile "$fa" \
        --organism euk \
        --output_dir "$outdir" \
        --format none \
        --mode slow-sequential

    gff="$outdir/output.gff3"

    if [ ! -f "$gff" ]; then
        echo "  ADVERTENCIA: No se encontró output.gff3 para $fa, se omite paso de corte."
        continue
    fi

    ########################################
    # 2) Extraer coordenadas de péptido señal del GFF3
    ########################################
    tmp_table=$(mktemp)

    # Filtra solo las líneas que corresponden a "signal_peptide"
    # Columnas: seqid (col1), tipo de feature (col3), inicio (col4), fin (col5)
    awk -F'\t' '
        $0 !~ /^#/ && $3 == "signal_peptide" {
            print $1 "\t" $4 "\t" $5
        }
    ' "$gff" > "$tmp_table"

    # Si no se detectan péptidos señal, saltar este archivo
    if [ ! -s "$tmp_table" ]; then
        echo "  No se detectaron péptidos señal en $fa."
        rm -f "$tmp_table"
        continue
    fi

    ########################################
    # 3) Crear FASTA cleavado + resumen
    ########################################
    cleaved_fasta="${base_noext}_signalp.fasta"  # FASTA con secuencias cleavadas
    summary_file="${core}_signalP_summary.tsv"   # Resumen de coordenadas

    # Código Python embebido para procesar el FASTA y crear salida
    python <<PY
import textwrap

fasta_path = "$fa"
table_path = "$tmp_table"
out_fasta_path = "$cleaved_fasta"
summary_path = "$summary_file"

# Cargar posiciones de corte
cleavage = {}
with open(table_path) as t:
    for line in t:
        line = line.strip()
        if not line:
            continue
        seqid, sp_start, sp_end = line.split("\t")
        if seqid not in cleavage:
            cleavage[seqid] = (int(sp_start), int(sp_end))

# Función para leer FASTA
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
                header_full = line[1:]
                header_id = header_full.split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header_id is not None:
        yield header_id, header_full, "".join(seq_lines)

# Abrir archivos de salida
with open(out_fasta_path, "w") as out_fa, open(summary_path, "w") as out_sum:
    out_sum.write("seqid\tsp_start\tsp_end\tcleavage_after_aa\toriginal_len\tnew_len\n")

    for sid, full_header, seq in read_fasta(fasta_path):
        if sid not in cleavage:
            continue

        sp_start, sp_end = cleavage[sid]
        cut_pos = sp_end  # índice 1-based
        new_seq = seq[cut_pos:]  # cortar SP

        if len(new_seq) == 0:
            continue

        # Escribir secuencia cleavada
        out_fa.write(f">{full_header} | signalp_cleaved\n")
        for chunk in textwrap.wrap(new_seq, 60):
            out_fa.write(chunk + "\n")

        # Escribir línea de resumen
        out_sum.write(f"{sid}\t{sp_start}\t{sp_end}\t{cut_pos}\t{len(seq)}\t{len(new_seq)}\n")
PY

    # Limpiar archivo temporal
    rm -f "$tmp_table"

    echo "  -> FASTA cleavado: $cleaved_fasta"
    echo "  -> Archivo resumen: $summary_file"
    echo "  -> Resultados crudos en: $outdir"
done

# Mensaje si no se procesó ningún archivo
if [ "$processed_any" = false ]; then
    echo "No se procesaron archivos FASTA. Verifique que filtered_outputs/ contenga .fasta/.fa/.faa"
    echo "o que los nombres proporcionados como argumentos existan en esa carpeta."
fi

echo "Proceso completado."
