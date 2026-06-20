#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import pandas as pd


# ============
# ARGUMENTOS
# ============

parser = argparse.ArgumentParser(
    description="Integración de taxonomía + genes VFDB/CARD + MGE"
)

parser.add_argument(
    "--sample",
    required=True,
    help="ID de la muestra"
)

parser.add_argument(
    "--taxonomy-report",
    required=True,
    help="Archivo kraken.report"
)

parser.add_argument(
    "--taxonomy-out",
    required=True,
    help="Archivo kraken.out"
)

parser.add_argument(
    "--vfdb",
    required=True,
    help="Resultados VFDB filtrados"
)

parser.add_argument(
    "--card",
    required=True,
    help="Resultados CARD filtrados"
)

parser.add_argument(
    "--mge",
    required=True,
    help="Resultados MGE filtrados"
)

parser.add_argument(
    "--output",
    required=True,
    help="Tabla final integrada"
)

args = parser.parse_args()


# ============
# FUNCIONES
# ============
def load_kraken_taxonomy(kraken_out):
    """
    Carga taxonomía por contig desde kraken.out
    """

    taxonomy_dict = {}

    with open(kraken_out) as file:

        for line in file:

            cols = line.strip().split("\t")

            # Validación básica
            if len(cols) < 3:
                continue

            contig = cols[1]
            taxid = cols[2]

            taxonomy_dict[contig] = taxid

    return taxonomy_dict

def load_kraken_report(kraken_report):
    """
    Carga nombres taxonómicos, ranks y dominio desde kraken.report
    """

    report_dict = {}

    current_domain = "Unknown"

    with open(kraken_report) as file:

        for line in file:

            cols = line.strip().split("\t")

            if len(cols) < 6:
                continue

            rank = cols[3].strip()
            taxid = cols[4].strip()
            name = cols[5].strip()

            # Detectar dominio
            if "Bacteria" in name:
                current_domain = "Bacteria"
                
            elif "Viruses" in name:
                current_domain = "Virus"
                
            elif "Fungi" in name:
                current_domain = "Fungi"
                
            elif "Archaea" in name:
                current_domain = "Archaea"

            report_dict[taxid] = {
                "taxonomy": name,
                "rank": rank,
                "domain": current_domain
            }
    
    return report_dict

def load_annotation(annotation_file, gene_type):
    """
    Carga VFDB o CARD y estandariza columnas
    """

    df = pd.read_csv(
        annotation_file,
        sep="\t"
    )

    # Si solo existe header y no hay hits
    if df.empty:
        return pd.DataFrame()

    # Selección de columans importantes
    df = df[[
        "#FILE",
        "SEQUENCE",
        "GENE",
        "START",
        "END",
        "%COVERAGE",
        "%IDENTITY",
        "DATABASE"
    ]]

    # Renombrar columnas
    df.columns = [
        "file",
        "contig",
        "element",
        "start",
        "end",
        "coverage",
        "identity",
        "database"
    ]

    # Añadir tipo
    df["type"] = gene_type

    return df

def load_mge(mge_file):
    """
    Carga resultados MobileElementFinder
    """

    df = pd.read_csv(
        mge_file,
        sep="\t",
    )

    if df.empty:
        return pd.DataFrame()

    df = df[[
        "contig",
        "name",
        "start",
        "end",
        "coverage",
        "identity",
        "type"
    ]]

    df.columns = [
        "contig",
        "element",
        "start",
        "end",
        "coverage",
        "identity",
        "mge_class"
    ]

    df["database"] = "MobileElementFinder"
    df["type"] = "MGE"

    return df


# ============
# CARGA DATOS
# ============

taxonomy_dict = load_kraken_taxonomy(args.taxonomy_out)

report_dict = load_kraken_report(args.taxonomy_report)

vfdb_df = load_annotation(
    args.vfdb,
    "virulence"
)

card_df = load_annotation(
    args.card,
    "AMR"
)

mge_df = load_mge(
    args.mge
)


# ============
# MARCAR ASOCIACIÓN CON MGE
# ============

vfdb_df["mge_associated"] = False
card_df["mge_associated"] = False

for _, mge in mge_df.iterrows():

    contig = mge["contig"]

    mge_start = mge["start"]
    mge_end = mge["end"]

    mask_vfdb = (
        (vfdb_df["contig"] == contig) &
        (vfdb_df["start"] <= mge_end) &
        (vfdb_df["end"] >= mge_start)
    )

    vfdb_df.loc[
        mask_vfdb,
        "mge_associated"
    ] = True

    mask_card = (
        (card_df["contig"] == contig) &
        (card_df["start"] <= mge_end) &
        (card_df["end"] >= mge_start)
    )

    card_df.loc[
        mask_card,
        "mge_associated"
    ] = True

# ============
# INTEGRACIÓN
# ============

vfdb_df["mge_class"] = pd.NA
card_df["mge_class"] = pd.NA
mge_df["mge_associated"] = False

final_df = pd.concat(
    [vfdb_df, card_df, mge_df],
    ignore_index=True
)

# Añadir muestra
final_df["sample"] = args.sample

# taxid por contig
final_df["taxid"] = final_df["contig"].map(
    lambda x: taxonomy_dict.get(x, "NA")
)

# Nombre taxonómico
final_df["taxonomy"] = final_df["taxid"].map(
    lambda x: report_dict.get(x, {}).get("taxonomy", "Unclassified")
)

# Rank taxonómico
final_df["rank"] = final_df["taxid"].map(
    lambda x: report_dict.get(x, {}).get("rank", "U")
)

# Dominio taxonómico
final_df["domain"] = final_df["taxid"].map(
    lambda x: report_dict.get(x, {}).get("domain", "Unknown")
)


# ============
# ORDEN COLUMNAS
# ============

final_df = final_df[[
    "sample",
    "contig",
    "domain",
    "taxonomy",
    "rank",
    "element",
    "type",
    "mge_class",
    "mge_associated",
    "coverage",
    "identity",
    "database"
]]


# ============
# OUTPUT
# ============

output_dir = os.path.dirname(args.output)

os.makedirs(
    output_dir,
    exist_ok=True
)

final_df.to_csv(
    args.output,
    sep="\t",
    index=False
)

print(f"[OK] Tabla integrada generada: {args.output}")
