# Pipeline bioinformático para la evaluación del potencial patogénico mediante análisis metagenómico

---

## 1.Descripción general

Este proyecto desarrolla un pipeline bioinformático modular y reproducible para el análisis de muestras metagenómicas bacterianas a partir de datos de secuenciación paired-end en formato FASTQ.

El objetivo principal es detectar e integrar información taxonómica y funcional relacionada con:

- Genes de virulencia.
- Genes asociados a resistencia antimicrobiana (AMR).
- Elementos genéticos móviles (MGEs).

A partir de esta información, el pipeline calcula métricas heurísticas de potencial patogénico que permite priorizar taxones y muestras según la presencia de determinantes funcionales potencialmente relevantes.

El proyecto forma parte de un Trabajo Fin de Máster (TFM) en Bioinformática.

---

## 2.Flujos de trabajo

El pipeline implementa las siguientes etapas:

FASTQ 

│ 

├── Control de calidad (FastQC) 

├── Filtrado y trimming (fastp) 

├── Ensamblado metagenómico (MEGAHIT) 

├── Clasificación taxonómica (Kraken2) 

├── Detección de genes de virulencia (VFDB) 

├── Detección de genes AMR (CARD) 

├── Detección de MGEs (MobileElementFinder) 

├── Filtrado de anotaciones 

├── Integración taxonómica-funcional 

├── Sistema de scoring heurístico 

└── Generación automática de reportes

---

## 3.Características principales

- Pipeline modular basado en Bash y Python.
- Configuración centralizada mediante archivos YAML.
- Ejecución reproducible mediante entornos Conda.
- Clasificaciónt taxonómica de contigs ensamblados.
- Detección de genes de virulencia y AMR.
- Detección de MGEs.
- Asociación automática entre genes funcionales y MGEs mediante solapamiento de coordenadas.
- Sistema de scoring heurístico configurable.
- Generación automática de tablas resumen y gráficos.
- Registro completo de ejecuciones mediante manifest y runlog.

---

## 4.Datos de entrada

**Origen:**

- Datos metagenómicos: Seqence Read Archive (SRA)

**Directorio de entrada:**

`data-raw`

**Formato aceptado**

Lecturas metagenómicas paired-end:

`sample_R1.fastq.gz`
`sample_R2.fastq.gz`

Cada muestra debe contener ambos archivos correspondientes a las lecturas forward y reverse.

**Restricciones/ética:**

- Datos de acceso público
- No contienen información sensible
- Los datos originales no se modifican

---

## 5.Herramientas utilizadas

| Etapa | Herramienta |
| :--- | :--- |
| Control de calidad | FastQC |
| Filtrado | fastp |
| Resumen QC | MultiQC |
| Ensamblado | MEGAHIT |
| Clasificación taxonómica | Kraken2 |
| Genes de virulencia | ABRicate + VFDB |
| Genes MAR | ABRicate + CARD |
| Detección de MGEs | MobileElementFinder |
| Integración y scoring | Python + Pandas |
| Reporting | Python + Pandas + Matplotlib |

---

## 6.Bases de datos utilizadas

**VFDB**

Base de datos especializada en factores de virulencia bacterianos.

**CARD**

Base de datos especializada en genes asociados a AMR.

**Kraken2 Standard Database**

Base de datos taxonómica utilizada para clasificación de contigs.

**MobileElementFinder Database**

Base de datos utilizada para la detección de elementos genéticos móviles.

---

## 7.Sistema de scoring

El pipeline incorpora un sistema heurístico diseñado para integrar la información funcional detectada.

Para cada gen identificado se calcula una puntuación basada en:

- Porcentaje de identidad.
- Porcentaje de cobertura.
- Tipo funcional.
- Asociación con MGEs.

Cada muestra y taxón son clasificados automáticamente en categorías:

- LOW
- MEDIUM
- HIGH

Los puesos y umbrales empleados pueden modificarse desde el archivo de configuración.

---

## 8.Salidas principales

**Resultados taxonómicos**

`taxonomy/`

- Clasificaciones Kraken2.
- Resúmenes taxonómicos.

**Resultados funcionales**

`annotation/`

- Genes de virulencia.
- Genes AMR.
- MGEs detectados.

**Integración**

`integration/`

- Tabla integrada final.
- Asociación genes-MGE.

**Scoring**

`scoring/`

- Scores por taxón.
- Scores globales.

**Reportes**

`report/`

- Summary.txt
- Top taxa
- Top genes
- Top MGEs
- Gráficos automatizados

---

## 9.Configuración

La configuración del pipeline se centraliza en:

`config/default.yaml`

Entre los parámetros configurables se incluyen:

- Rutas de bases de datos.
- Enotornos Conda.
- Umbrales de identidad.
- Umbrales de cobertura.
- Pesos del sistema de scoring.
- Parámetros de visualización.
- Opciones de ejecución.

---

## 10.Reproducibilidad

El pipeline incorpora distintos mecanismos destinados a garantizar la reproducibilidad:

- Configuración centralizada mediante YAML.
- Entornos Conda específicos para cada módulo.
- Registro de ejecuciones mediante runlog.
- Control de muestras mediante manifest.
- Conservación de resultados intermedios.
- Posibilidad de reanudar ejecuciones desde etapas concretas.

---

## 11.Requisitos

- Linux
- Conda o Miniconda
- Python 3
- Bases de datos previamente instaladas

---

## 12.Cómo ejecutar

**Comando de ejemplo**

```bash
bash scripts/metagenomic_pipeline.sh -i data-raw/META001/ -o results/exp01
```

---

## 13.Estado del proyecto

Versión funcional completa.

Actualmente implementa:

- Análisis taxonómico.
- Anotación funcional.
- Detección de MGEs.
- Integración taxonómica-funcional.
- Scoring heurístico.
- Generación automática de reportes.
