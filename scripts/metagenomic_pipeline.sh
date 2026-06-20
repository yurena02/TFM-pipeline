#!/usr/bin/env bash

export LC_ALL=C

help(){
	cat << EOF

Uso:
	$0 -i one_metagenomes_dir -o output_dir [opciones]
	$0 -I multiples_metagenomes_dir -o output_dir [opciones]

Descripción:
Pipeline metagenómico para la detección de genes de virulencia y resistencia antimicrobiana con integración en un sistema de scoring.

Argumentos obligatorios:
	-o DIR Directorio de salida

Entrada (elegir una):
	-i DIR Carpeta de una muestra (ej: META001/)
	-I DIR Carpeta con múltiples muestras

Opciones:
	-c FILE Archivo de configuración YAML (default: config/default.yaml)
	-s STEP Ejecutar desde este paso (qc, assembly, taxonomy, mapping, annotation, filtering, integration, scoring, report)
	-e STEP Ejecutar solo este paso
	-u STEP Ejecutar hasta este paso
	-t INT Número de threads (default: 4)
	-h Mostrar esta ayuda

Formato esperado por muestra:
Cada muestra debe contener archivos FASTQ paired-end:
	sample_R1.fastq
	sample_R2.fastq

Ejemplos:
	Ejecutar pipeline completo con una muestra:
	$0 -i data-raw/META001/ -o results/

	Ejecutar múltiples muestras:
	$0 -I data-raw/ -o results/

	Ejecutar desde anotación:
	$0 -i data-raw/META001/ -o results/ -s annotation

	Ejecutar solo scoring:
	$0 -i data-raw/META001/ -o results/ -e scoring

EOF
exit
}

# ======== DEFAULTS  ========

COMMAND="$(printf '%q ' "$0" "$@")"
PIPELINE_STATUS="OK"
FAILED_STEP="none"
SAMPLE_DIR=""
INPUT_DIR=""
OUTPUT=""
CONFIG=""
START_STEP=""
ONLY_STEP=""
END_STEP=""
THREADS=""
SAMPLES=()
RUN_PIPELINE=false
CURRENT_STEP="initialization"
YAML_ENV="yaml_env"

# ======== GETOPTS ========

while getopts "i:I:o:c:s:e:u:t:h" option; do
	case $option in
		i) SAMPLE_DIR=$OPTARG;;
		I) INPUT_DIR=$OPTARG;;
		o) OUTPUT=$OPTARG;;
		c) CONFIG=$OPTARG;;
		s) START_STEP=$OPTARG;;
		e) ONLY_STEP=$OPTARG;;
		u) END_STEP=$OPTARG;;
		t) THREADS=$OPTARG;;
		h) help;;
	esac
done


# ======== VALIDACION DE ERRORES ========

set -Eeuo pipefail

trap 'PIPELINE_STATUS="ERROR"; FAILED_STEP="$CURRENT_STEP"' ERR

# Input obligatorio (una de las dos opciones)
if [[ -z "$SAMPLE_DIR" && -z "$INPUT_DIR" ]]; then
	echo "Error; debe proporcionar -i o -I"
	help; exit 1
fi

# No permitir ambos input
if [[ -n "$SAMPLE_DIR" && -n "$INPUT_DIR" ]]; then
	echo "Error: no se pueden usar -i y -I simultáneamente"
	exit 1
fi

# Existencia del input
if [[ -n "$SAMPLE_DIR" ]]; then
	[[ ! -d "$SAMPLE_DIR" ]] && echo "Error: no existe $SAMPLE_DIR" && exit 1
fi

if [[ -n "$INPUT_DIR" ]]; then
	[[ ! -d "$INPUT_DIR" ]] && echo "Error: no existe $INPUT_DIR" && exit 1
fi

# Existencia de R1 y R2 para -i

if [[ -n "$SAMPLE_DIR" ]]; then
	R1=$(ls "$SAMPLE_DIR"/*R1*.fastq* 2>/dev/null || true)
	R2=$(ls "$SAMPLE_DIR"/*R2*.fastq* 2>/dev/null || true)

	[[ -z "$R1" || -z "$R2" ]] && {
		echo "Error: no se encontraron archivos R1/R2 en $SAMPLE_DIR"
		exit 1
	}
fi

# Output obligatorio
[[ -z "$OUTPUT" ]] && echo "Error: debe proporcionar -o" && help && exit 1

mkdir -p "$OUTPUT"

# Valores por defecto
[[ -z "$CONFIG" ]] && CONFIG="config/default.yaml"
[[ -z "$THREADS" ]] && THREADS=4

# Validación de que threads sea un número
[[ "$THREADS" =~ ^[0-9]+$ ]] || {
	echo "Error: threads debe ser un número entero"
	exit 1
}

# No permitir combinaciones conflictivas
if [[ -n "$START_STEP" && -n "$ONLY_STEP" ]]; then
	echo "Error: no se pueden usar -s y -e simultáneamente"
	exit 1
fi

if [[ -n "$START_STEP" && -n "$END_STEP" ]]; then
	echo "Error: no se pueden usar -u y -s simultáneamente"
	exit 1
fi

if [[ -n "$ONLY_STEP" && -n "$END_STEP" ]]; then
	echo "Error: no se pueden usar -u y -e simultáneamente"
	exit 1
fi

# Validar pasos
VALID_STEPS=(
	qc
	assembly
	taxonomy
	annotation
	filtering
	integration
	scoring
	report
)

check_step(){
	local step=$1
	for s in "${VALID_STEPS[@]}"; do
		[[ "$s" == "$step" ]] && return 0
	done
	echo "Error: paso inválido -> $step"
	exit 1
}

[[ -n "$START_STEP" ]] && check_step "$START_STEP"
[[ -n "$ONLY_STEP" ]] && check_step "$ONLY_STEP"
[[ -n "$END_STEP" ]] && check_step "$END_STEP"


# ======== RUTAS GENERALES ========

set_sample_paths(){

	local sample_path=$1

	SAMPLE=$(basename "$sample_path")

	SAMPLE_OUT="$OUTPUT/$SAMPLE"

	INTERMEDIATE_DIR="$SAMPLE_OUT/intermediate"
	FINAL_RESULTS_DIR="$SAMPLE_OUT/final_results"

	QC_DIR="$INTERMEDIATE_DIR/QC"
	ASSEMBLY_DIR="$INTERMEDIATE_DIR/assembly"
	TAXONOMY_DIR="$INTERMEDIATE_DIR/taxonomy"
	ANNOTATION_DIR="$INTERMEDIATE_DIR/annotation"
	
	VIRULENCE_DIR="$ANNOTATION_DIR/virulence"
	AMR_DIR="$ANNOTATION_DIR/amr"
	MGE_DIR="$ANNOTATION_DIR/mge"

	INTEGRATION_DIR="$FINAL_RESULTS_DIR/integration"
	SCORING_DIR="$FINAL_RESULTS_DIR/scoring"
	REPORT_DIR="$FINAL_RESULTS_DIR/report"
}

RUNLOG="$OUTPUT/runlog.txt"
MANIFEST="docs/manifest-resultados.tsv"


# ======== HELPERS ========

validate_file(){

	local file=$1
	local message=$2

	if [[ ! -f "$file" ]]; then
		echo "Error: $message"
		exit 1
	fi
}

validate_nonempty(){

	local file=$1
	local message=$2

	if [[ ! -s "$file" ]]; then
		echo "Error: $message"
		exit 1
	fi
}

detect_reads(){
	
	local dir=$1

	R1=$(ls "$dir"/*R1*.fastq* 2>/dev/null || true)
	R2=$(ls "$dir"/*R2*.fastq* 2>/dev/null || true)
}

get_config(){

	local key=$1

	conda run -n "$YAML_ENV" \
		python scripts/load_config.py \
		"$CONFIG" \
		"$key"
}

log_step(){
	echo ""
	echo "==> Paso $1 iniciado"
}

log_sample(){
	echo ""
	echo "Procesando muestra: $SAMPLE"
}

log_fin_step(){
	echo ""
	echo "==> $1 completado"
}

# ======== CONFIG.YAML ========

QC_TOOLS=$(get_config conda_envs.qc)
FASTP_QUALITY=$(get_config qc.fastp.qualified_quality_phred)
FASTP_LENGTH=$(get_config qc.fastp.length_required)
FASTP_ADAPTER=$(get_config qc.fastp.detect_adapter_for_pe)
MULTIQC_ENV=$(get_config conda_envs.multiqc)

ASSEMBLY_ENV=$(get_config conda_envs.assembly)
MEGAHIT_MIN_CONTIG=$(get_config assembly.megahit.min_contig_len)
MEGAHIT_PRESET=$(get_config assembly.megahit.presets)

KRAKEN_DB=$(get_config databases.kraken2.db)
TAXONOMY_ENV=$(get_config conda_envs.taxonomy)
KRAKEN_MEMORY_MAPPING=$(get_config taxonomy.kraken2.memory_mapping)

VFDB_DB=$(get_config databases.abricate.vfdb)
CARD_DB=$(get_config databases.abricate.card)
MGE_DB=$(get_config databases.mobile_element_finder.db)

ANNOTATION_ENV=$(get_config conda_envs.annotation)
MGE_ENV=$( get_config conda_envs.mge)

VFDB_MIN_IDENTITY=$(get_config filtering.vfdb.min_identity)
VFDB_MIN_COVERAGE=$(get_config filtering.vfdb.min_coverage)
CARD_MIN_IDENTITY=$(get_config filtering.card.min_identity)
CARD_MIN_COVERAGE=$(get_config filtering.card.min_coverage)
MGE_MIN_IDENTITY=$(get_config filtering.mge.min_identity)
MGE_MIN_COVERAGE=$(get_config filtering.mge.min_coverage)

FILTERING_ENV=$(get_config conda_envs.filtering)

INTEGRATION_ENV=$(get_config conda_envs.integration)

SCORING_ENV=$(get_config conda_envs.scoring)
VIRULENCE_WEIGHT=$(get_config scoring.weights.virulence)
AMR_WEIGHT=$(get_config scoring.weights.amr)
MGE_WEIGHT=$(get_config scoring.weights.mge)
LOW_THRESHOLD=$(get_config scoring.thresholds.low)
HIGH_THRESHOLD=$(get_config scoring.thresholds.high)

REPORT_ENV=$(get_config conda_envs.report)
TOP_TAXA=$(get_config report.top_taxa)
TOP_GENES=$(get_config report.top_genes)
TOP_MGE=$(get_config report.top_mge)
PLOT_DPI=$(get_config report.plots.dpi)
BARPLOT_WIDTH=$(get_config report.plots.barplot_width)
BARPLOT_HEIGHT=$(get_config report.plots.barplot_height)
PIECHART_WIDTH=$(get_config report.plots.piechart_width)
PIECHART_HEIGHT=$(get_config report.plots.piechart_height)

# ======== INPUT (FASTQ) ========

if [[ -n "$SAMPLE_DIR" ]]; then
        SAMPLES+=("$SAMPLE_DIR")
fi

if [[ -n "$INPUT_DIR" ]]; then
        for dir in "$INPUT_DIR"/*/; do
                [[ -d "$dir" ]] && SAMPLES+=("$dir")
        done
fi

# Validar que haya muestras

[[ ${#SAMPLES[@]} -eq 0 ]] && {
        echo "Error: no se encontraron muestras"
        exit 1
}

# Existencia de R1 y R2

VALID_SAMPLES=()

for SAMPLE_PATH in "${SAMPLES[@]}"; do
        detect_reads "$SAMPLE_PATH"

        if [[ -z "$R1" || -z "$R2" ]]; then
                echo "Error: no se encontraron R1/R2 en $SAMPLE_PATH"
                echo "Muestra inválida, se omite"
                continue
        fi

        VALID_SAMPLES+=("$SAMPLE_PATH")

done

# Reemplazar array original por muestras válidas
SAMPLES=("${VALID_SAMPLES[@]}")

# Verificar que quede al menos una muestra válida
[[ ${#SAMPLES[@]} -eq 0 ]] && {
        echo "Error: no hay muestras válidas"
        exit 1
}


# ======== MANIFEST HELPER ========
mkdir -p docs

RUN_ID=$(basename "$OUTPUT")

# Registrar outputs generados en el manifest global
if [[ ! -f "$MANIFEST" ]]; then
        echo -e "run_id\tsample_id\tstep\tsubstep\tfile_type\tfile_path\tfile_name\tcreated_at\ttool\tdescription" > "$MANIFEST"

fi

log_manifest(){

        local sample_id=$1
        local step=$2
        local substep=$3
        local file_type=$4
        local file_path=$5
        local file_name=$6
        local tool=$7
        local description=$8

        local created_at=$(date "+%Y-%m-%d %H:%M:%S")

        echo -e "${RUN_ID}\t${sample_id}\t${step}\t${substep}\t${file_type}\t${file_path}\t${file_name}\t${created_at}\t${tool}\t${description}" >> "$MANIFEST"
}


# ======== RUNLOG HELPER ========

write_runlog(){

	[[ -z "${OUTPUT:-}" ]] && return

        END_TIME=$(date +%s)
        END_HUMAN=$(date)

        DURATION=$((END_TIME - START_TIME))

        MINUTES=$((DURATION / 60))
        SECONDS_REMAINING=$((DURATION % 60))

        {
                echo "ID de ejecución: $RUN_ID"
                echo ""

                echo "Fecha inicio: $START_HUMAN"
                echo "Fecha fin: $END_HUMAN"
                echo "Duración: ${MINUTES}m${SECONDS_REMAINING}s"
                echo ""

                echo "Entrada:"
                [[ -n "$SAMPLE_DIR" ]] && echo "$SAMPLE_DIR"
                [[ -n "$INPUT_DIR" ]] && echo "$INPUT_DIR"
                echo ""

                echo "Salida:"
                echo "$OUTPUT"
                echo ""

                echo "Comando:"
                echo "bash $COMMAND"
                echo ""

                echo "Parámetros:"
                printf "\t- threads: %s\n" "$THREADS"
                printf "\t- config: %s\n" "$CONFIG"

                [[ -n "$START_STEP" ]] && \
                        printf "\t- start_step: %s\n" "$START_STEP"

                [[ -n "$ONLY_STEP" ]] && \
                        printf "\t- only_step: %s\n" "$ONLY_STEP"

                [[ -n "$END_STEP" ]] && \
                        printf "\t- end_step: %s\n" "$END_STEP"

                echo ""

                echo "Versiones herramientas:"

                printf "\t- fastqc: %s\n" \
                        "$(conda run -n "$QC_TOOLS" fastqc --version 2>/dev/null | head -n1 || echo 'N/A')"
                
                printf "\t- fastp: %s\n" \
                        "$(conda run -n "$QC_TOOLS" fastp --version 2>/dev/null | head -n1 || echo 'N/A')"

                printf "\t- multiqc: %s\n" \
                        "$(conda run -n "$MULTIQC_ENV" multiqc --version 2>/dev/null | head -n1 || echo 'N/A')"

                printf "\t- megahit: %s\n" \
                        "$(conda run -n "$ASSEMBLY_ENV" megahit --version 2>/dev/null | head -n1 || echo 'N/A')"

                printf "\t- kraken2: %s\n" \
                        "$(conda run -n "$TAXONOMY_ENV" kraken2 --version 2>/dev/null | head -n1 || echo 'N/A')"

                printf "\t- abricate: %s\n" \
                        "$(conda run -n "$ANNOTATION_ENV" abricate --version 2>/dev/null | head -n1 || echo 'N/A')"

		printf "\t- mobileelementfinder: %s\n" \
			"$(conda run -n "$MGE_ENV" mefinder --version 2>/dev/null | head -n1 || echo 'N/A')"

                printf "\t- python: %s\n" \
                        "$(conda run -n "$INTEGRATION_ENV" python --version 2>/dev/null | head -n1 || echo 'N/A')"

                printf "\t- pandas: %s\n" \
                        "$(conda run -n "$INTEGRATION_ENV" python -c 'import pandas as pd; print(pd.__version__)' 2>/dev/null || echo 'N/A')"

                printf "\t- matplotlib: %s\n" \
                        "$(conda run -n "$REPORT_ENV" python -c 'import matplotlib; print(matplotlib.__version__)' 2>/dev/null || echo 'N/A')"

                printf "\t- pyyaml: %s\n" \
                        "$(conda run -n "$YAML_ENV" python -c 'import yaml; print(yaml.__version__)' 2>/dev/null || echo 'N/A')"
                echo ""

                echo "Entorno:"
                printf "\t- OS: %s\n" "$(lsb_release -ds)"
                echo ""

                echo "Versión del código:"
                printf "\t- commit: %s\n" "$(git rev-parse --short HEAD 2>/dev/null || echo 'N/A')"
                printf "\t- branch: %s\n" "$(git branch --show-current 2>/dev/null || echo 'N/A')"
                echo ""

                echo "Resultado: $PIPELINE_STATUS"

                if [[ "$PIPELINE_STATUS" == "ERROR" ]]; then
                        echo "Paso fallido: $FAILED_STEP"
                fi
		echo ""
		echo ""

        } >> "$RUNLOG"
}


# ======== CHECK TOOLS ========

check_tools(){
	echo "==> Comprobando herramientas"
	
	local missing=false
	
	declare -A TOOLS=(
		["$QC_TOOLS"]="fastqc fastp"
		["$MULTIQC_ENV"]="multiqc"
		["$ASSEMBLY_ENV"]="megahit"
		["$TAXONOMY_ENV"]="kraken2"
		["$ANNOTATION_ENV"]="abricate"
		["$MGE_ENV"]="mefinder blastn"
		["$FILTERING_ENV"]="abricate"
		["$INTEGRATION_ENV"]="python"
		["$SCORING_ENV"]="python"
		["$REPORT_ENV"]="python"
		["$YAML_ENV"]="python"
	)
	
	declare -A PYTHON_LIBS=(
		["$INTEGRATION_ENV"]="pandas"
		["$SCORING_ENV"]="pandas"
		["$REPORT_ENV"]="pandas matplotlib"
		["$YAML_ENV"]="yaml"
	)
	
	for env in "${!TOOLS[@]}"; do
		
		# Verificar que existe el entorno
		if ! conda env list | grep -q "^$env "; then
			echo "Error: entorno Conda no encontrado -> $env"
			missing=true
			continue
		fi
		
		# Verificar herramientas dentro del entorno
		for tool in ${TOOLS[$env]}; do
			
			if ! conda run -n "$env" which "$tool" >/dev/null 2>&1; then
				echo "Error: herramienta no encontrada -> $tool ($env)"
				missing=true
			else
				echo "OK: $tool ($env)"
			fi
		done
	done
	
	for env in "${!PYTHON_LIBS[@]}"; do
		
		# Verificar librerías de Python
		for lib in ${PYTHON_LIBS[$env]}; do
			if ! conda run -n "$env" python -c "import $lib" >/dev/null 2>&1; then
				echo "Error: librería no encontrada -> $lib ($env)"
				missing=true
			else
				echo "OK: $lib ($env)"
			fi
		done
	done
	
	if [[ "$missing" == true ]]; then
		echo "==> Faltan herramientas, librerías o entornos"
		exit 1
	fi
	
	echo "==> Todas las herramientas disponibles"
}


# ======== M1. CONTROL DE CALIDAD ========

run_qc(){
	log_step "QC"

	for SAMPLE_PATH in "${SAMPLES[@]}"; do

		set_sample_paths "$SAMPLE_PATH"
		log_sample

		# Definir rutas
		FASTQC_RAW="$QC_DIR/fastqc_raw"
		FASTQC_TRIMMED="$QC_DIR/fastqc_trimmed"
		FASTP_OUT="$QC_DIR/fastp"

		# Crear las carpetas de salida si no existen
		mkdir -p "$FASTQC_RAW" "$FASTQC_TRIMMED" "$FASTP_OUT"


		# Detectar FASTQ
		detect_reads "$SAMPLE_PATH"
		
		[[ -z "$R1" || -z "$R2" ]] && {
			echo "Error: no se encontraron R1/R2 en $SAMPLE_PATH"
			continue
		}

		# ==== FASTQC RAW ====
		echo "FastQC (raw)"

		conda run -n "$QC_TOOLS" fastqc -t "$THREADS" "$R1" "$R2" -o "$FASTQC_RAW"

		# Manifest
		log_manifest \
			"$SAMPLE" \
			"qc" \
			"fastqc_raw" \
			"html" \
			"$FASTQC_RAW/" \
			"$(basename "$R1" .fastq.gz)_fastqc.html" \
			"fastqc" \
			"FastQC report (raw R1)"

		log_manifest \
                        "$SAMPLE" \
                        "qc" \
                        "fastqc_raw" \
                        "html" \
                        "$FASTQC_RAW/" \
                        "$(basename "$R2" .fastq.gz)_fastqc.html" \
                        "fastqc" \
                        "FastQC report (raw R2)"

		# ==== FASTP ====

		echo "Trimming con fastp"

		BASE_R1=$(basename "$R1" .fastq.gz)
		BASE_R2=$(basename "$R2" .fastq.gz)

		OUT_R1="$FASTP_OUT/${BASE_R1}_trimmed.fastq.gz"
		OUT_R2="$FASTP_OUT/${BASE_R2}_trimmed.fastq.gz"

		FASTP_EXTRA=""

		[[ "$FASTP_ADAPTER" == "True" ]] && FASTP_EXTRA="--detect_adapter_for_pe"
		conda run -n "$QC_TOOLS" fastp \
			-i "$R1" -I "$R2" \
			-o "$OUT_R1" -O "$OUT_R2" \
			-h "$FASTP_OUT/${SAMPLE}.html" \
			-j "$FASTP_OUT/${SAMPLE}.json" \
			-w "$THREADS" \
			--qualified_quality_phred "$FASTP_QUALITY" \
			--length_required "$FASTP_LENGTH" \
			$FASTP_EXTRA

		# Manifest
		log_manifest \
			"$SAMPLE" \
			"qc" \
			"fastp" \
			"fastq.gz" \
			"$FASTP_OUT/" \
			"$(basename "$OUT_R1")" \
			"fastp" \
			"Trimmed reads R1"

		log_manifest \
                        "$SAMPLE" \
                        "qc" \
                        "fastp" \
                        "fastq.gz" \
                        "$FASTP_OUT/" \
                        "$(basename "$OUT_R2")" \
                        "fastp" \
                        "Trimmed reads R2"

		log_manifest \
                        "$SAMPLE" \
                        "qc" \
                        "fastp" \
                        "json" \
                        "$FASTP_OUT/" \
			"${SAMPLE}.json" \
                        "fastp" \
                        "Fastp stats JSON"

		log_manifest \
                        "$SAMPLE" \
                        "qc" \
                        "fastp" \
                        "html" \
                        "$FASTP_OUT/" \
                        "${SAMPLE}.html" \
                        "fastp" \
                        "Fastp report HTML"

		# ==== FASTQC TRIMMED ====
		echo "FastQC (trimmed)"
		conda run -n "$QC_TOOLS" fastqc -t "$THREADS" "$OUT_R1" "$OUT_R2" -o "$FASTQC_TRIMMED"

		# Manifest
		log_manifest \
                        "$SAMPLE" \
                        "qc" \
                        "fastqc_trimmed" \
                        "html" \
                        "$FASTQC_TRIMMED/" \
			"$(basename "$OUT_R1" .fastq.gz)_fastqc.html" \
                        "fastqc" \
                        "FastQC report (trimmed R1)"

		log_manifest \
                        "$SAMPLE" \
                        "qc" \
                        "fastqc_trimmed" \
                        "html" \
                        "$FASTQC_TRIMMED/" \
			"$(basename "$OUT_R2" .fastq.gz)_fastqc.html" \
                        "fastqc" \
                        "FastQC report (trimmed R2)"

	done

	# ==== MULTIQC GLOBAL ====
	MULTIQC_OUT="$OUTPUT/multiqc"
	mkdir -p "$MULTIQC_OUT"
	
	echo "Ejecutando MultiQC"

	conda run -n "$MULTIQC_ENV" multiqc "$OUTPUT" -o "$MULTIQC_OUT"

	# Manifest
	TOTAL_SAMPLES=${#SAMPLES[@]}

        log_manifest \
                        "GLOBAL" \
                        "qc" \
                        "multiqc" \
                        "html" \
                        "$MULTIQC_OUT/" \
                        "multiqc_report.html" \
                        "multiqc" \
                        "Aggregated QC report (${TOTAL_SAMPLES})"

	log_fin_step "QC"
}



# ======== M2. ASSEMBLY -> CONTIGS ========

run_assembly(){

	log_step "ASSEMBLY"
	
	for SAMPLE_PATH in "${SAMPLES[@]}"; do

		set_sample_paths "$SAMPLE_PATH"
		log_sample

		# Definir rutas
		FASTP_OUT="$QC_DIR/fastp"
		MEGAHIT_OUT="$ASSEMBLY_DIR/final.contigs.fa"
		CONTIGS="$ASSEMBLY_DIR/contigs.fa"

		# Skip si assembly ya existe
		if [[ -s "$CONTIGS" ]]; then
			echo "Assembly ya existe para $SAMPLE -> saltando"
			continue
		fi

		# Si existe la carpeta pero no el resultado -> limpiarla
		if [[ -d "$ASSEMBLY_DIR" ]]; then
			echo "Directorio assembly existe pero incompleto -> limpiando"
			rm -rf "$ASSEMBLY_DIR"
		fi

		# Detectar reads filtrados
		if [[ ! -d "$QC_DIR" ]]; then
			echo "Error: QC no ejecutado para $SAMPLE"
			exit 1
		fi

		R1=$(ls "$FASTP_OUT"/*R1*trimmed.fastq.gz 2>/dev/null || true)
		R2=$(ls "$FASTP_OUT"/*R2*trimmed.fastq.gz 2>/dev/null || true)

		if [[ -z "$R1" || -z "$R2" ]]; then
			echo "Error: no se encontraron reads filtrados en $QC_DIR"
			continue
		fi

		# Ejecución de MEGAHIT
		echo "Ejecutando MEGAHIT"

		MEGAHIT_EXTRA=""

		[[ -n "$MEGAHIT_PRESET" ]] && MEGAHIT_EXTRA="--presets $MEGAHIT_PRESET"

		conda run -n "$ASSEMBLY_ENV" megahit \
			-1 "$R1" \
			-2 "$R2" \
			-o "$ASSEMBLY_DIR" \
			-t "$THREADS" \
			--min-contig-len "$MEGAHIT_MIN_CONTIG" \
			$MEGAHIT_EXTRA

		validate_file "$MEGAHIT_OUT" "MEGAHIT no generó contigs para $SAMPLE"

		# Renombrar para estandarizar
		mv "$MEGAHIT_OUT" "$CONTIGS"

		# Manifest
		log_manifest \
			"$SAMPLE" \
			"assembly" \
			"megahit_contigs" \
			"fasta" \
			"$ASSEMBLY_DIR/" \
			"contigs.fa" \
			"megahit" \
			"Final assembled contigs"

		echo "Assembly completado: $ASSEMBLY_DIR/contigs.fa"
	done

	log_fin_step "ASSEMBLY"
}


# ======== M3. CLASIFICACIÓN TAXONÓMICA (CONTIGS) ========

run_taxonomy(){

	log_step "TAXONOMY"

	# Validar base de datos
	if [[ ! -d "$KRAKEN_DB" ]]; then
		echo "Error: no se encontró la base de datos de Kraken2 en $KRAKEN_DB"
		exit 1
	fi

	for SAMPLE_PATH in "${SAMPLES[@]}"; do

		set_sample_paths "$SAMPLE_PATH"
		log_sample

		# Definir rutas
		CONTIGS="$ASSEMBLY_DIR/contigs.fa"
		KRAKEN_OUT="$TAXONOMY_DIR/kraken.out"
		KRAKEN_REPORT="$TAXONOMY_DIR/kraken.report"

		mkdir -p "$TAXONOMY_DIR"

		# Skip si ya existe
		if [[ -s "$KRAKEN_REPORT" ]]; then
			echo "Taxonomía ya existe para $SAMPLE -> saltando"
			continue
		fi

		# Validar input
		validate_file "$CONTIGS" "no se encontró contigs.fa para $SAMPLE"

		echo "Ejecutando Kraken 2"

		KRAKEN_EXTRA=""
		
		[[ "$KRAKEN_MEMORY_MAPPING" == "True" ]] && KRAKEN_EXTRA="--memory-mapping"

		conda run -n "$TAXONOMY_ENV" kraken2 \
			--db "$KRAKEN_DB" \
			--threads "$THREADS" \
			--report "$KRAKEN_REPORT" \
			--output "$KRAKEN_OUT" \
			$KRAKEN_EXTRA \
			"$CONTIGS"

		# Validación output
		validate_nonempty "$KRAKEN_REPORT" "Kraken2 no generó resultados para $SAMPLE"

		# Manifest
		log_manifest \
			"$SAMPLE" \
			"taxonomy" \
			"kraken2_report" \
			"report" \
			"$TAXONOMY_DIR/" \
			"kraken.report" \
			"kraken2" \
			"Taxonomic classification report"

		echo "Taxonomía completada: $KRAKEN_REPORT"

	done

	log_fin_step "TAXONOMY"
}

# ======== M4. ANOTACIONES FUNCIONALES ========

# ==== VIRULENCE ====

run_virulence(){

	echo "==> Anotación virulencia iniciada"

	# Base de datos VFDB para ABRicate

	for SAMPLE_PATH in "${SAMPLES[@]}"; do

		set_sample_paths "$SAMPLE_PATH"
		log_sample

		# Rutas
		CONTIGS="$ASSEMBLY_DIR/contigs.fa"
		VFDB_RAW="$VIRULENCE_DIR/vfdb_raw.tsv"
		
		mkdir -p "$VIRULENCE_DIR"

		# Skip si ya existe el archivo
		if [[ -s "$VFDB_RAW" ]]; then
			echo "Anotación ya existe para $SAMPLE -> saltando"
			continue
		fi

		# Validar input
		validate_file "$CONTIGS" "no se encontró contigs.fa para $SAMPLE"

		# == ABRICATE RAW ==
		
		echo "Ejecutando ABRicate (VFDB)"

		conda run -n "$ANNOTATION_ENV" abricate \
			--db "$VFDB_DB" \
			"$CONTIGS" \
			> "$VFDB_RAW"

		# Validar output
		validate_nonempty "$VFDB_RAW" "ABRicate no generó resultados para $SAMPLE"

		# Manifest
		log_manifest \
			"$SAMPLE" \
			"annotation" \
			"vfdb_raw" \
			"tsv" \
			"$VIRULENCE_DIR/" \
			"vfdb_raw.tsv" \
			"abricate" \
			"Raw VFDB annotation results"	

		 echo "Anotación virulencia completada: $SAMPLE"

	done

}

# ==== AMR ====

run_amr(){

	echo "==> Anotación AMR iniciada"

	for SAMPLE_PATH in "${SAMPLES[@]}"; do

		set_sample_paths "$SAMPLE_PATH"

		# Rutas
		CONTIGS="$ASSEMBLY_DIR/contigs.fa"
		CARD_RAW="$AMR_DIR/card_raw.tsv"
		
		mkdir -p "$AMR_DIR"

		# Skip si ya existe el archivo
		if [[ -s "$CARD_RAW" ]]; then
			echo "AMR ya existe para $SAMPLE -> saltando"
			continue
		fi

		# Validar input
		validate_file "$CONTIGS" "no se encontró contigs.fa para $SAMPLE"

		# == ABRICATE RAW ==
		echo "Ejecutando ABRicate (CARD)"

		conda run -n "$ANNOTATION_ENV" abricate \
			--db "$CARD_DB" \
			"$CONTIGS" \
			> "$CARD_RAW"

		# Validar output
		validate_nonempty "$CARD_RAW" "ABRicate no generó resultados AMR para $SAMPLE"

		# Manifest
		log_manifest \
			"$SAMPLE" \
			"annotation" \
			"card_raw" \
			"tsv" \
			"$AMR_DIR/" \
			"card_raw.tsv" \
			"abricate" \
			"Raw CARD annotation results"

		echo "AMR completado: $SAMPLE"

	done
}

# ==== MOBILE GENETIC ELEMENTS (MGE) ====

run_mge(){

	# Validar base de datos
        if [[ ! -d "$MGE_DB" ]]; then
                echo "Error: no se encontró la base de datos de MGE en $MGE_DB"
                exit 1
        fi

	echo "==> Anotación MGE iniciada"

		for SAMPLE_PATH in "${SAMPLES[@]}"; do

		set_sample_paths "$SAMPLE_PATH"
		log_sample

		# Rutas
		CONTIGS="$ASSEMBLY_DIR/contigs.fa"

		MGE_PREFIX="$MGE_DIR/mobile_element_finder"
		MGE_CSV="$MGE_DIR/mobile_element_finder.csv"
		MGE_RAW="$MGE_DIR/mge_raw.tsv"

		MGE_TMP="/tmp/mge_finder"

		mkdir -p "$MGE_DIR"

		# Skip si ya existe
		if [[ -s "$MGE_RAW" ]]; then
			echo "MGE ya existe para $SAMPLE -> saltando"
			continue
		fi

		# Validar input
		validate_file "$CONTIGS" "no se encontró contigs.fa para $SAMPLE"

		# Borrar caché si existe
		if [ -d "$MGE_TMP" ]; then
			echo "[INFO] Limpiando caché temporal de MobileElementFinder"
			rm -rf "$MGE_TMP"
		fi

		echo "Ejecutando MobileElementFinder"

		conda run -n "$MGE_ENV" mefinder find \
			-c "$CONTIGS" \
			--db-path "$MGE_DB" \
			-t "$THREADS" \
			"$MGE_PREFIX"
		
		validate_nonempty "$MGE_CSV" "MobileElementFinder no generó resultados para $SAMPLE"

		tr ',' '\t' < "$MGE_CSV" | grep -v '^#' > "$MGE_RAW"
		
		validate_nonempty "$MGE_RAW" "No se ha creado $MGE_RAW para $SAMPLE"

		# Manifest
		log_manifest \
			"$SAMPLE" \
			"annotation" \
			"mge_raw" \
			"tsv" \
			"$MGE_DIR/" \
			"mge_raw.tsv" \
			"MobileElementFinder" \
			"Raw mobile genetic element predictions"

		echo "MGE completado: $SAMPLE"

	done
}


# ==== EJECUCIÓN COMPLETA DE ANOTACIONES FUNCIONALES ====

run_annotation(){

	log_step "ANNOTATION"

	run_virulence
	run_amr
	run_mge

	log_fin_step "ANNOTATION"
}

# ======== M5. FILTRADO (THRESHOLDS + LIMPIEZA) ========

run_filtering(){

	log_step "FILTERING"

	# Thresholds globales

	for SAMPLE_PATH in "${SAMPLES[@]}"; do

		set_sample_paths "$SAMPLE_PATH"
		log_sample

		# Rutas
		VFDB_RAW="$VIRULENCE_DIR/vfdb_raw.tsv"
		CARD_RAW="$AMR_DIR/card_raw.tsv"
		MGE_RAW="$MGE_DIR/mge_raw.tsv"

		VFDB_FILTERED="$VIRULENCE_DIR/vfdb_filtered.tsv"
		CARD_FILTERED="$AMR_DIR/card_filtered.tsv"
		MGE_FILTERED="$MGE_DIR/mge_filtered.tsv"

		VFDB_SUMMARY="$VIRULENCE_DIR/vfdb_summary.tsv"
		CARD_SUMMARY="$AMR_DIR/card_summary.tsv"
		MGE_SUMMARY="$MGE_DIR/mge_summary.tsv"

		# Validar input
		validate_file "$VFDB_RAW" "No se encontró vfdb_raw.tsv para $SAMPLE"
		validate_file "$CARD_RAW" "No se encontró card_raw.tsv para $SAMPLE"
		validate_file "$MGE_RAW" "No se encontró mge_raw.tsv para $SAMPLE"
		# ==== VFDB ====

		echo "Filtrando VFDB"

		awk -F '\t' '
			BEGIN {OFS="\t"}
			NR==1 {print; next}
			($10 >= '$VFDB_MIN_COVERAGE' && $11 >= '$VFDB_MIN_IDENTITY')
                ' "$VFDB_RAW" > "$VFDB_FILTERED"

		# Summary VFDB 

		echo "Generando resumen VFDB"

                if [[ $(wc -l < "$VFDB_FILTERED") -gt 1 ]]; then
			conda run -n "$FILTERING_ENV" \
			       abricate --summary "$VFDB_FILTERED" > "$VFDB_SUMMARY"
                else
			cp "$VFDB_FILTERED" "$VFDB_SUMMARY"
                fi

                # Manifest

		log_manifest \
                        "$SAMPLE" \
                        "filtering" \
                        "vfdb_filtered" \
                        "tsv" \
                        "$VIRULENCE_DIR/" \
                        "vfdb_filtered.tsv" \
                        "awk" \
                        "Filtered VFDB results (identity >= ${VFDB_MIN_IDENTITY} and coverage >= ${VFDB_MIN_COVERAGE})"

                log_manifest \
			"$SAMPLE" \
			"filtering" \
                        "vfdb_summary" \
                        "tsv" \
                        "$VIRULENCE_DIR/" \
                        "vfdb_summary.tsv" \
                        "abricate" \
                        "VFDB summary table"

		# ==== CARD ====

		echo "Filtrando CARD"

		awk -F '\t' '
			BEGIN {OFS="\t"}
                        NR==1 {print; next}
                        ($10 >= '$CARD_MIN_COVERAGE' && $11 >= '$CARD_MIN_IDENTITY')
		' "$CARD_RAW" > "$CARD_FILTERED"

		# Summary CARD

		echo "Generando resumen CARD"

                if [[ $(wc -l < "$CARD_FILTERED") -gt 1 ]]; then
			conda run -n "$FILTERING_ENV" \
				abricate --summary "$CARD_FILTERED" > "$CARD_SUMMARY"
                else
			cp "$CARD_FILTERED" "$CARD_SUMMARY"
                fi

                # Manifest
                
		log_manifest \
                        "$SAMPLE" \
                        "filtering" \
                        "card_filtered" \
                        "tsv" \
                        "$AMR_DIR/" \
                        "card_filtered.tsv" \
                        "awk" \
                        "Filtered CARD results (identity >= ${CARD_MIN_IDENTITY} and coverage >= ${CARD_MIN_COVERAGE})"

		log_manifest \
                        "$SAMPLE" \
                        "filtering" \
                        "card_summary" \
                        "tsv" \
                        "$AMR_DIR/" \
                        "card_summary.tsv" \
                        "abricate" \
                        "CARD summary table"

		# ==== MGE ====

		echo "Filtrando MGE"

		LC_ALL=C awk -F '\t' \
			-v min_id="$MGE_MIN_IDENTITY" \
			-v min_cov="$MGE_MIN_COVERAGE" '
			BEGIN {OFS="\t"}
			NR==1 {print; next}
			($9+0 >= min_id+0 && $10+0 >= min_cov+0) {
				split($13,a," ")
				$13=a[1]
				print
			}
		' "$MGE_RAW" > "$MGE_FILTERED"

		# Summary MGE

		echo "Generando resumen MGE"

		awk -F '\t' '
			BEGIN {OFS="\t"}
			!/^#/ && $1!="mge_no" {
				count[$5]++
			}
			END {
				print "type","count"
				for (t in count)
					print t,count[t]
			}
			' "$MGE_FILTERED" > "$MGE_SUMMARY"

		# Manifest

		log_manifest \
			"$SAMPLE" \
			"filtering" \
			"mge_filtered" \
			"tsv" \
			"$MGE_DIR/" \
			"mge_filtered.tsv" \
			"awk" \
			"Filtered MobileElementFinder results"

		log_manifest \
			"$SAMPLE" \
			"filtering" \
			"mge_summary" \
			"tsv" \
			"$MGE_DIR/" \
			"mge_summary.tsv" \
			"awk" \
			"Summary of predicted mobile genetic elements by type"

		echo "Filtering completado: $SAMPLE"

	done

	log_fin_step "FILTERING"
}

# ======== M6. INTEGRACIÓN DE RESULTADOS ========

run_integration(){

	log_step "INTEGRATION"

	for SAMPLE_PATH in "${SAMPLES[@]}"; do

		set_sample_paths "$SAMPLE_PATH"
		log_sample

                # Rutas
                VFDB_FILTERED="$VIRULENCE_DIR/vfdb_filtered.tsv"
                CARD_FILTERED="$AMR_DIR/card_filtered.tsv"
		MGE_FILTERED="$MGE_DIR/mge_filtered.tsv"

		KRAKEN_OUT="$TAXONOMY_DIR/kraken.out"
		KRAKEN_REPORT="$TAXONOMY_DIR/kraken.report"

		INTEGRATION_OUT="$INTEGRATION_DIR/final_table.tsv"

		mkdir -p "$INTEGRATION_DIR"

		# Skip si ya existe
		if [[ -s "$INTEGRATION_OUT" ]]; then
			echo "Integration ya existe para $SAMPLE -> saltando"
			continue
		fi

		# Validar input
		
		validate_file "$VFDB_FILTERED" "No se encontró vfdb_filtered.tsv para $SAMPLE"
		validate_file "$CARD_FILTERED" "No se encontró card_filtered.tsv para $SAMPLE"
		validate_file "$MGE_FILTERED" "No se encontró mge_filtered.tsv para $SAMPLE"
		validate_file "$KRAKEN_OUT" "No se encontró kraken.out para $SAMPLE"
		validate_file "$KRAKEN_REPORT" "No se encontró kraken.report para $SAMPLE"

		# Ejecución
		
		echo "Ejecutando integración"

		conda run -n "$INTEGRATION_ENV" python scripts/integration.py \
			--sample "$SAMPLE" \
			--taxonomy-report "$KRAKEN_REPORT" \
			--taxonomy-out "$KRAKEN_OUT" \
			--vfdb "$VFDB_FILTERED" \
			--card "$CARD_FILTERED" \
			--mge "$MGE_FILTERED" \
			--output "$INTEGRATION_OUT"

		# Validar output
		validate_nonempty "$INTEGRATION_OUT" "No se generó final_table.tsv para $SAMPLE"

		# Manifest
		log_manifest \
			"$SAMPLE" \
			"integration" \
			"final_table" \
			"tsv" \
			"$INTEGRATION_DIR/" \
			"final_table.tsv" \
			"pandas" \
			"Integrated taxonomy + virulence + AMR + MGE table"

		echo "Integración completada: $SAMPLE"

	done

	log_fin_step "INTEGRATION"
}

# ======== M7. SCORING ========

run_scoring(){

	log_step "SCORING"

	# Output global compartido
	GLOBAL_RESULTS_DIR="$OUTPUT/final_results"
	mkdir -p "$GLOBAL_RESULTS_DIR"

	GLOBAL_SCORING="$GLOBAL_RESULTS_DIR/global_scoring.tsv"

	for SAMPLE_PATH in "${SAMPLES[@]}"; do

		set_sample_paths "$SAMPLE_PATH"
		log_sample

		# Rutas
		INTEGRATED_TABLE="$INTEGRATION_DIR/final_table.tsv"
		TAXONOMY_SCORING="$SCORING_DIR/taxonomy_scoring.tsv"

		mkdir -p "$SCORING_DIR"

		# Skip si ya existe
		if [[ -s "$TAXONOMY_SCORING" && -s "$GLOBAL_SCORING" ]]; then
			echo "Scoring ya existe para $SAMPLE -> saltando"
			continue
		fi

		# Validar input
		validate_file "$INTEGRATED_TABLE" "No se encontró final_table.tsv para $SAMPLE"

		# Ejecución
		echo "Ejecutando scoring"

		conda run -n "$SCORING_ENV" python scripts/scoring.py \
			--sample "$SAMPLE" \
			--integrated-table "$INTEGRATED_TABLE" \
			--sample-output "$TAXONOMY_SCORING" \
			--global-output "$GLOBAL_SCORING" \
			--virulence-weight "$VIRULENCE_WEIGHT" \
			--amr-weight "$AMR_WEIGHT" \
			--mge-weight "$MGE_WEIGHT" \
			--low-threshold "$LOW_THRESHOLD" \
			--high-threshold "$HIGH_THRESHOLD"

		# Validar outputs
		validate_nonempty "$TAXONOMY_SCORING" "No se generó taxonomy_scoring.tsv para $SAMPLE"

		validate_nonempty "$GLOBAL_SCORING" "No se generó global_scoring.tsv"

		# Manifest
		log_manifest \
			"$SAMPLE" \
			"scoring" \
			"taxonomy_scoring" \
			"tsv" \
			"$SCORING_DIR/" \
			"taxonomy_scoring.tsv" \
			"pandas" \
			"Taxonomic pathogenicity scoring table"

		log_manifest \
			"$SAMPLE" \
			"scoring" \
			"global_scoring" \
			"tsv" \
			"$OUTPUT/" \
			"global_scoring.tsv" \
			"pandas" \
			"Global pathogenicity scoring summary"

		echo "Scoring completado: $SAMPLE"

	done

	log_fin_step "SCORING"
}


# ======== M8. REPORT ========

run_report(){

	log_step "REPORT"

	for SAMPLE_PATH in "${SAMPLES[@]}"; do

		set_sample_paths "$SAMPLE_PATH"
		log_sample

		# Rutas
		INTEGRATED_TABLE="$INTEGRATION_DIR/final_table.tsv"
		GLOBAL_SCORING="$OUTPUT/final_results/global_scoring.tsv"
		TAXONOMY_SCORING="$SCORING_DIR/taxonomy_scoring.tsv"

		mkdir -p "$REPORT_DIR"

		# Outputs
		SUMMARY_FILE="$REPORT_DIR/summary.txt"
		TOP_TAXA_FILE="$REPORT_DIR/top_taxa.tsv"
		TOP_GENES_FILE="$REPORT_DIR/top_genes.tsv"
		TOP_MGES_FILE="$REPORT_DIR/top_mges.tsv"
		TAXA_BARPLOT_FILE="$REPORT_DIR/taxa_barplot.png"
		FUNCTIONAL_PIECHART_FILE="$REPORT_DIR/functional_piechart.png"
		MGE_BARPLOT_FILE="$REPORT_DIR/mge_barplot.png"

		# Skip si ya existen outputs
		if [[ -s "$SUMMARY_FILE" && -s "$TOP_TAXA_FILE" && -s "$TOP_GENES_FILE" && -s "$TOP_MGES_FILE" && -s "$TAXA_BARPLOT_FILE" && -s "$FUNCTIONAL_PIECHART_FILE" && -s "$MGE_BARPLOT_FILE" ]]; then
			echo "Report ya existe para $SAMPLE -> saltando"
			continue
		fi

		# Validar inputs
		validate_file "$INTEGRATED_TABLE" "No se encontró final_table.tsv para $SAMPLE"
		validate_file "$GLOBAL_SCORING" "No se encontró global_scoring.tsv para $SAMPLE"
		validate_file "$TAXONOMY_SCORING" "No se encontró taxonomy_scoring.tsv para $SAMPLE"

		# Ejecución
		echo "Generando report"

		conda run -n "$REPORT_ENV" python scripts/report.py \
			--sample "$SAMPLE" \
			--integrated-table "$INTEGRATED_TABLE" \
			--global-scoring "$GLOBAL_SCORING" \
			--taxonomy-scoring "$TAXONOMY_SCORING" \
			--output-dir "$REPORT_DIR" \
			--top-taxa "$TOP_TAXA" \
			--top-genes "$TOP_GENES" \
			--top-mge "$TOP_MGE" \
			--plot-dpi "$PLOT_DPI" \
			--barplot-width "$BARPLOT_WIDTH" \
			--barplot-height "$BARPLOT_HEIGHT" \
			--piechart-width "$PIECHART_WIDTH" \
			--piechart-height "$PIECHART_HEIGHT"

		# Validar outputs
		validate_nonempty "$SUMMARY_FILE" "No se generó summary.txt para $SAMPLE"
		validate_nonempty "$TOP_TAXA_FILE" "No se generó top_taxa.tsv para $SAMPLE"
		validate_nonempty "$TOP_GENES_FILE" "No se generó top_genes.tsv para $SAMPLE"
		validate_nonempty "$TOP_MGES_FILE" "No se generó top_mges.tsv para $SAMPLE"
		validate_nonempty "$TAXA_BARPLOT_FILE" "No se generó taxa_barplot.png para $SAMPLE"
		validate_nonempty "$FUNCTIONAL_PIECHART_FILE" "No se generó functional_piechart.png para $SAMPLE"
		validate_nonempty "$MGE_BARPLOT_FILE" "No se generó mge_barplot.png para $SAMPLE"

		# Manifest
		log_manifest \
			"$SAMPLE" \
			"report" \
			"summary" \
			"txt" \
			"$REPORT_DIR/" \
			"summary.txt" \
			"pandas" \
			"General report summary"

		log_manifest \
			"$SAMPLE" \
			"report" \
			"top_taxa" \
			"tsv" \
			"$REPORT_DIR/" \
			"top_taxa.tsv" \
			"pandas" \
			"Top taxa ranked by pathogenicity score"

		log_manifest \
			"$SAMPLE" \
			"report" \
			"top_genes" \
			"tsv" \
			"$REPORT_DIR/" \
			"top_genes.tsv" \
			"pandas" \
			"Top detected virulence and AMR genes"

		log_manifest \
			"$SAMPLE" \
			"report" \
			"top_mges" \
			"tsv" \
			"$REPORT_DIR/" \
			"top_mges.tsv" \
			"pandas" \
			"Top detected MGE elements"

		log_manifest \
			"$SAMPLE" \
			"report" \
			"taxa_barplot" \
			"png" \
			"$REPORT_DIR/" \
			"taxa_barplot.png" \
			"matplotlib" \
			"Barplot of top taxa scores"

		log_manifest \
			"$SAMPLE" \
			"report" \
			"functional_piechart" \
			"png" \
			"$REPORT_DIR/" \
			"functional_piechart.png" \
			"matplotlib" \
			"Pie chart of risk categories"

		log_manifest \
			"$SAMPLE" \
			"report" \
			"mge_barplot" \
			"png" \
			"$REPORT_DIR/" \
			"mge_barplot.png" \
			"matplotlib" \
			"Barplot of top MGE elements"

		echo "Reporte completado: $SAMPLE"

	done

	log_fin_step "REPORT"
}


# ======== EJECUCIONES DEL SCRIPT ========

# Separar la ejecución completa de la ejecución por partes (-s, -e y -u)
run_step(){

	local step_name=$1
	local step_function=$2

	# echo "DEBUG: ONLY_STEP='$ONLY_STEP' STEP='$step_name'"
	
	CURRENT_STEP="$step_name"

	# Si se ha definido ONLY_STEP -> solo ejecutar ese paso
	if [[ -n "$ONLY_STEP" ]]; then
		if [[ "$step_name" == "$ONLY_STEP" ]]; then
			echo "==> Ejecutando solo el paso $step_name"
			$step_function
		fi
		return
	fi

	# Si se ha definido START_STEP -> activar ejecución desde ahí
	if [[ -n "$START_STEP" ]]; then
		if [[ "$step_name" == "$START_STEP" ]]; then
			RUN_PIPELINE=true
		fi

		if [[ "${RUN_PIPELINE}" == true ]]; then
			echo "==> Ejecutando paso: $step_name"
			$step_function
		fi
		return
	fi

	# Si se ha definido END_STEP -> ejecutar hasta ese paso
	if [[ -n "$END_STEP" ]]; then
		echo "==> Ejecutando paso: $step_name"
		$step_function

		if [[ "$step_name" == "$END_STEP" ]]; then
			echo "==> END_STEP alcanzado ($END_STEP), deteniendo pipeline"
			exit 0
		fi
		return
	fi

	# Caso normal -> ejecutar todo
	echo "==> Ejecutando paso: $step_name"
	$step_function
}


# Flujo de ejecución completo
run_pipeline(){
	
	echo "==> INICIO PIPELINE"

	run_step qc run_qc
	run_step assembly run_assembly
	run_step taxonomy run_taxonomy
	run_step annotation run_annotation
	run_step filtering run_filtering
	run_step integration run_integration
	run_step scoring run_scoring
	run_step report run_report

	echo "==> FIN PIPELINE"
}


# Ejecución del script
main(){
	START_TIME=$(date +%s)
	START_HUMAN=$(date)

	check_tools
	run_pipeline
}

trap write_runlog EXIT
main
