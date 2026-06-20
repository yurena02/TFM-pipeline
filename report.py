#!/usr/bin/env python
# coding: utf-8


import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


# ============
# ARGUMENTOS
# ============

parser = argparse.ArgumentParser(
    description="Generación de reportes finales"
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
    "--global-scoring",
    required=True,
    help="Tabla global_scoring.tsv"
)

parser.add_argument(
    "--taxonomy-scoring",
    required=True,
    help="Tabla taxonomy_scoring.tsv"
)

parser.add_argument(
    "--output-dir",
    required=True,
    help="Directorio de salida report/"
)

parser.add_argument(
    "--top-taxa",
    type=int,
    required=True,
)

parser.add_argument(
    "--top-genes",
    type=int,
    required=True,
)

parser.add_argument(
    "--top-mge",
    type=int,
    required=True
)

parser.add_argument(
    "--plot-dpi",
    type=int,
    required=True,
)

parser.add_argument(
    "--barplot-width",
    type=float,
    required=True,
)

parser.add_argument(
    "--barplot-height",
    type=float,
    required=True,
)

parser.add_argument(
    "--piechart-width",
    type=float,
    required=True,
)

parser.add_argument(
    "--piechart-height",
    type=float,
    required=True
)

args = parser.parse_args()

# ============
# PARÁMETROS
# ============

TOP_TAXA = args.top_taxa
TOP_GENES = args.top_genes
TOP_MGE = args.top_mge

PLOT_DPI = args.plot_dpi

BARPLOT_WIDTH = args.barplot_width
BARPLOT_HEIGHT = args.barplot_height

PIECHART_WIDTH = args.piechart_width
PIECHART_HEIGHT = args.piechart_height

# ============
# FUNCIONES
# ============

def normalize_taxonomy(taxonomy):

    if pd.isna(taxonomy):
        return "Unclassified"

    taxonomy = str(taxonomy).strip()

    unknown_values = [
        "",
        "NA",
        "Unknown",
        "Unclassified"
    ]

    if taxonomy in unknown_values:
        return "Unclassified"

    return taxonomy


def create_summary(
    sample,
    global_row,
    top_taxon,
    mge_total,
    mge_associated,
    output_file
):
    with open(output_file, "w") as f:

        f.write(
            "===============\n"
        )

        f.write(
            "METAGENOMIC PATHOGENICITY REPORT\n"
        )

        f.write(
            "===============\n\n"
        )

        f.write(f"Sample: {sample}\n\n")

        f.write("GLOBAL RESULTS\n")
        f.write("---------------\n")

        f.write(
            f"Virulence genes: "
            f"{global_row['virulence_genes']}\n"
        )

        f.write(
            f"AMR genes: "
            f"{global_row['amr_genes']}\n"
        )

        f.write(
            f"Detected MGEs: "
            f"{mge_total}\n"
        )

        f.write(
            f"MGE-associated elements: "
            f"{mge_associated}\n"
        )

        f.write(
            f"Virulence score: "
            f"{global_row['virulence_score']}\n"
        )

        f.write(
            f"AMR score: "
            f"{global_row['amr_score']}\n"
        )

        f.write(
            f"Total score: "
            f"{global_row['total_score']}\n"
        )

        f.write(
            f"Risk category: "
            f"{global_row['category']}\n\n"
        )

        f.write("TOP TAXON\n")
        f.write("-------------\n")

        f.write(
            f"Taxonomy: "
            f"{top_taxon['taxonomy']}\n"
        )

        f.write(
            f"Total score: "
            f"{top_taxon['total_score']}\n"
        )

        f.write(
            f"Virulence genes: "
            f"{top_taxon['virulence_genes']}\n"
        )

        f.write(
            f"AMR genes: "
            f"{top_taxon['amr_genes']}\n"
        )


# ============
# CARGA DE DATOS
# ============

integrated_df = pd.read_csv(
    args.integrated_table,
    sep="\t"
)

taxonomy_df = pd.read_csv(
    args.taxonomy_scoring,
    sep="\t"
)

global_df = pd.read_csv(
    args.global_scoring,
    sep="\t"
)

# Validación básica
if integrated_df.empty:
    raise ValueError(
        "La tabla integrada está vacía"
    )

if taxonomy_df.empty:
    raise ValueError(
        "La tabla taxonomy_scoring está vacía"
    )

if global_df.empty:
    raise ValueError(
        "La tabla global_scoring está vacía"
    )


# ============
# LIMPIEZA
# ============

integrated_df["taxonomy"] = integrated_df[
    "taxonomy"
].apply(normalize_taxonomy)

# Eliminar duplicados
integrated_df = integrated_df.drop_duplicates(
    subset=[
        "contig",
        "element",
        "type"
    ]
)



# ============
# CREAR DIRECTORIO
# ============

os.makedirs(
    args.output_dir,
    exist_ok=True
)


# ============
# TOP TAXA
# ============

top_taxa_df = taxonomy_df.sort_values(
    by="total_score",
    ascending=False
).head(TOP_TAXA)

top_taxa_out = os.path.join(
    args.output_dir,
    "top_taxa.tsv"
)

top_taxa_df.to_csv(
    top_taxa_out,
    sep="\t",
    index=False
)


# ============
# TOP GENES
# ============

top_genes_df = integrated_df.groupby(
    ["element", "type"]
).size().reset_index(name="count")

top_genes_df = top_genes_df.sort_values(
    by="count",
    ascending=False
).head(TOP_GENES)

top_genes_out = os.path.join(
    args.output_dir,
    "top_genes.tsv"
)

top_genes_df.to_csv(
    top_genes_out,
    sep="\t",
    index=False
)


# ============
# TOP MGES
# ============

top_mges_df = (
    integrated_df[
        integrated_df["mge_class"].notna()
    ]
    .groupby("mge_class")
    .size()
    .reset_index(name="count")
    .sort_values(
        by="count",
        ascending=False
    )
    .head(TOP_MGE)
)

top_mges_out = os.path.join(
    args.output_dir,
    "top_mges.tsv"
)

top_mges_df.to_csv(
    top_mges_out,
    sep="\t",
    index=False
)


# ============
# SUMARY TXT
# ============

global_row = global_df[
    global_df["sample"] == args.sample
].iloc[0]

top_taxon = top_taxa_df.iloc[0]

summary_out = os.path.join(
    args.output_dir,
    "summary.txt"
)

mge_total = len(
    integrated_df[
        integrated_df["type"] == "MGE"
    ]
)

mge_associated = len(
    integrated_df[
        integrated_df["mge_associated"] == True
    ]
)

create_summary(
    args.sample,
    global_row,
    top_taxon,
    mge_total,
    mge_associated,
    summary_out
)


# ============
# BARPLOT TAXONÓMICO
# ============

barplot_out = os.path.join(
    args.output_dir,
    "taxa_barplot.png"
)

plot_df = top_taxa_df.head(TOP_TAXA)

plt.figure(figsize=(BARPLOT_WIDTH, BARPLOT_HEIGHT))

plt.bar(
    plot_df["taxonomy"],
    plot_df["total_score"]
)

plt.xticks(
    rotation=45,
    ha="right"
)

plt.ylabel("Total score")

plt.xlabel("Taxonomy")

plt.title(
    f"Top taxa by pathogenicity score"
    f"({args.sample})"
)

plt.tight_layout()

plt.savefig(
    barplot_out,
    dpi=PLOT_DPI
)

plt.close()


# ============
# PIECHART FUNCIONAL
# ============

piechart_out = os.path.join(
    args.output_dir,
    "functional_piechart.png"
)

virulence_count = len(
    integrated_df[
        integrated_df["type"] == "virulence"
    ]
)

amr_count = len(
    integrated_df[
        integrated_df["type"] == "AMR"
    ]
)


plt.figure(figsize=(PIECHART_WIDTH, PIECHART_HEIGHT))

plt.pie(
    [virulence_count, amr_count],
    labels=["Virulence", "AMR"],
    autopct="%1.1f%%"
)

plt.title(
    f"Functional composition ({args.sample})"
)
plt.tight_layout()

plt.savefig(
    piechart_out,
    dpi=PLOT_DPI
)

plt.close()


# ============
# MGE BARPLOT
# ============

mge_barplot_out = os.path.join(
    args.output_dir,
    "mge_barplot.png"
)

mge_plot_df = integrated_df[
    integrated_df["type"] == "MGE"
]

if not mge_plot_df.empty:
    
    mge_plot_df = (
        mge_plot_df
        .groupby("mge_class")
        .size()
        .reset_index(name="count")
        .sort_values(
            by="count",
            ascending=False
        )
    )

    plt.figure(
        figsize=(
            BARPLOT_WIDTH,
            BARPLOT_HEIGHT
        )
    )

    plt.bar(
        mge_plot_df["mge_class"],
        mge_plot_df["count"]
    )

    plt.xticks(
        rotation=45,
        ha="right"
    )

    plt.ylabel(
        "Count"
    )

    plt.xlabel(
        "MGE class"
    )

    plt.title(
        f"MGE composition ({args.sample})"
    )

    plt.tight_layout()

    plt.savefig(
        mge_barplot_out,
        dpi=PLOT_DPI
    )

    plt.close()


# ============
# OUTPUT FINAL
# ============

print(
    f"[OK] Reporte generado para {args.sample}"
)

print(
    f"[OK] Summary: {summary_out}"
)

print(
    f"[OK] Top taxa: {top_taxa_out}"
)

print(
    f"[OK] Top genes: {top_genes_out}"
)

print(
    f"[OK] Top MGEs: {top_mges_out}"
)

print(
    f"[OK] Barplot: {barplot_out}"
)

print(
    f"[OK] Piechart: {piechart_out}"
)

print(
    f"[OK] MGE barplot: {mge_barplot_out}"
)
