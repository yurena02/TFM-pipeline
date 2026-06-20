#!/usr/bin/env python
# coding: utf-8


import os
import argparse
import pandas as pd


# ============
# ARGUMENTOS
# ============

parser = argparse.ArgumentParser(
    description="Scoring de potencial patogénico"
)

parser.add_argument(
    "--sample",
    required=True,
    help="ID de la muestra"
)

parser.add_argument(
    "--integrated-table",
    required=True,
    help="Tabla integrada final_table.tsv"
)

parser.add_argument(
    "--sample-output",
    required=True,
    help="Output scoring por taxonomía"
)

parser.add_argument(
    "--global-output",
    required=True,
    help="Output scoring global"
)

parser.add_argument(
    "--virulence-weight",
    type=float,
    required=True
)

parser.add_argument(
    "--amr-weight",
    type=float,
    required=True
)

parser.add_argument(
    "--mge-weight",
    type=float,
    required=True
)

parser.add_argument(
    "--low-threshold",
    type=int,
    required=True
)

parser.add_argument(
    "--high-threshold",
    type=int,
    required=True
)

args = parser.parse_args()

# ============
# PARÁMETROS
# ============

# Pesos globales
VIRULENCE_WEIGHT = args.virulence_weight
AMR_WEIGHT = args.amr_weight
MGE_WEIGHT = args.mge_weight

# Categorías heurísticas
LOW_THRESHOLD = args.low_threshold
HIGH_THRESHOLD = args.high_threshold


# ============
# FUNCIONES
# ============

def normalize_taxonomy(taxonomy):
    """
    Normaliza taxonomías desconocidas
    """

    if pd.isna(taxonomy):
        return "Unclassified"

    taxonomy = str(taxonomy).strip()

    unknown_values = [
        "NA",
        "Unknown",
        "Unclassified",
        ""
    ]

    if taxonomy in unknown_values:
        return "Unclassified"

    return taxonomy


def calculate_element_score(row):
    """
    Calcula score individual por gen

    Fórmula:
    (coverage * identity) / 10000
    """

    return (
        row["coverage"] *
        row["identity"]
    ) / 10000


def assign_category(score):
    """
    Asignación heurística de categorías
    """

    if score < LOW_THRESHOLD:
        return "LOW"

    elif score < HIGH_THRESHOLD:
        return "MEDIUM"

    else:
        return "HIGH"

def apply_weight(row):

    if row["type"] == "virulence":
        score = row["element_score"] * VIRULENCE_WEIGHT

    elif row["type"] == "AMR":
        score = row["element_score"] * AMR_WEIGHT

    else:
        return 0

    if (
        row["type"] != "MGE"
        and row["mge_associated"] == True
    ):
        score *= MGE_WEIGHT

    return score


# ============
# CARGA DATOS
# ============

df = pd.read_csv(
    args.integrated_table,
    sep="\t"
)

# Validación básica
if df.empty:
    raise ValueError(
        f"No se encontraron datos en {args.integrated_table}"
    )


# ============
# LIMPIEZA
# ============

# Normalizar taxonomías
df["taxonomy"] = df["taxonomy"].apply(
    normalize_taxonomy
)

# Eliminar duplicados exactos
df = df.drop_duplicates(
    subset=[
        "contig",
        "element",
        "type"
    ]
)


# ============
# SCORE POR GEN
# ============

df["element_score"] = df.apply(
    calculate_element_score,
    axis=1
)

# Aplicar pesos
df["weighted_score"] = df.apply(
    apply_weight,
    axis=1
)

scoring_df = df[
    df["type"] != "MGE"
]

# ============
# SCORING POR TAXONOMÍA
# ============

species_rows = []

for taxonomy, group in scoring_df.groupby("taxonomy"):

    virulence_df = group[
        group["type"] == "virulence"
    ]

    amr_df = group[
        group["type"] == "AMR"
    ]

    mge_df = group[
        group["type"] == "MGE"
    ]

    virulence_genes = len(virulence_df)
    amr_genes = len(amr_df)

    virulence_score = virulence_df[
        "weighted_score"
    ].sum()

    amr_score = amr_df[
        "weighted_score"
    ].sum()

    total_score = (
        virulence_score +
        amr_score
    )

    species_rows.append({
        "sample": args.sample,
        "taxonomy": taxonomy,
        "virulence_genes": virulence_genes,
        "amr_genes": amr_genes,
        "virulence_score": round(
            virulence_score, 3
        ),
        "amr_score": round(
            amr_score, 3
        ),
        "total_score": round(
            total_score, 3
        ),
        "category": assign_category(
            total_score
        )
    })

species_df = pd.DataFrame(
    species_rows
)

species_df = species_df[
    species_df["total_score"] > 0
]

# Ordenar por score
species_df = species_df.sort_values(
    by="total_score",
    ascending=False
)


# ============
# SCORING GLOBAL
# ============
global_virulence_df = df[
    df["type"] == "virulence"
]

global_amr_df = df[
    df["type"] == "AMR"
]

global_virulence_genes = len(
    global_virulence_df
)

global_amr_genes = len(
    global_amr_df
)

global_virulence_score = (
    global_virulence_df[
        "weighted_score"
    ].sum()
)

global_amr_score = (
    global_amr_df[
        "weighted_score"
    ].sum()
)

global_total_score = (
    global_virulence_score +
    global_amr_score
)

global_df = pd.DataFrame([{
    "sample": args.sample,
    "virulence_genes": global_virulence_genes,
    "amr_genes": global_amr_genes,
    "virulence_score": round(
        global_virulence_score, 3
    ),
    "amr_score": round(
        global_amr_score, 3
    ),
    "total_score": round(
        global_total_score, 3
    ),
    "category": assign_category(
        global_total_score
    )
}])


# ============
# OUTPUTS
# ============

# Crear directorios
sample_output_dir = os.path.dirname(
    args.sample_output
)

global_output_dir = os.path.dirname(
    args.global_output
)

os.makedirs(
    sample_output_dir,
    exist_ok=True
)

os.makedirs(
    global_output_dir,
    exist_ok=True
)

# Exportar scoring por taxonomía
species_df.to_csv(
    args.sample_output,
    sep="\t",
    index=False
)

# Exportar scoring global
if os.path.exists(args.global_output):

    previous_df = pd.read_csv(
        args.global_output,
        sep="\t"
    )

    global_df = pd.concat(
        [previous_df, global_df],
        ignore_index=True
    )

# Evitar duplicados de muestra
global_df = global_df.drop_duplicates(
    subset=["sample"],
    keep="last"
)

global_df.to_csv(
    args.global_output,
    sep="\t",
    index=False
)

print(
    f"[OK] Scoring generado para {args.sample}"
)

print(
    f"[OK] Output taxonómico: "
    f"{args.sample_output}"
)

print(
    f"[OK] Output global: "
    f"{args.global_output}"
)
