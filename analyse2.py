"""Analysis of biochemical data from ABC efflux pump inhibitors.

This script retrieves and processes experimental data from PubChem
for three ABC transporter proteins (ABCB1, ABCG2, ABCC1) and
identifies their inhibitors through a multi-step cleaning pipeline.

Usage:
======
    python analyse.py

    The following CSV files must be present in the working directory:
        - ABCB1.csv
        - ABCG2.csv
        - ABCC1.csv

    Output files:
        - ABCB1_final.csv, ABCG2_final.csv, ABCC1_final.csv
        - distribution_activite.png
        - chevauchements.png
"""

__authors__ = "Ludine"
__contact__ = "viet-khoa.tran-nguyen@u-paris.fr"
__date__ = "2026"
__version__ = "1.0.0"

# Modules internes Python
import csv
from collections import Counter
from io import StringIO

# Modules externes
import matplotlib.pyplot as plt
import pandas as pd
import requests

# Constante globale : limite de taille des cellules CSV
# nécessaire car les citations dans PubChem sont très longues
CSV_FIELD_LIMIT = 10000000
csv.field_size_limit(CSV_FIELD_LIMIT)

# Constante : seuil d'activité en µM
ACTIVITY_THRESHOLD = 10.0

# Constante : nombre de colonnes attendues dans les fichiers PubChem
EXPECTED_COLUMNS = 32

# Constantes : AIDs pour la détection des faux positifs
AIDS_LUCIFERASE_COUNTER = [585, 485341]
AIDS_LUCIFERASE_CONFIRM = [584, 485294]
AIDS_AGGREGATORS = [411]
AIDS_AUTOFLUORESCENT = [587, 588, 590, 591, 592, 593, 594]

# Constante : colonnes à conserver dans les fichiers finaux
FINAL_COLUMNS = [
    "Substance_SID", "Compound_CID", "Activity",
    "Activity_Type", "Activity_Qualifier", "Activity_Value"
]


# ============================================================
# ÉTAPE 0 : CHARGEMENT DES FICHIERS CSV
# ============================================================

def charger_csv(nom_fichier):
    """Load a PubChem CSV file into a pandas DataFrame.

    PubChem CSV files have an unusual format: each data line is
    wrapped in double quotes and ends with semicolons. This function
    handles this format by cleaning each line before parsing.

    Parameters
    ----------
    nom_fichier : str
        Path to the CSV file to load.

    Returns
    -------
    pandas.DataFrame
        A cleaned DataFrame with 32 columns and stripped whitespace.
    """
    lignes_propres = []

    with open(nom_fichier, "r", encoding="utf-8") as f:
        for ligne in f:
            # Remove trailing newline and semicolons.
            ligne = ligne.rstrip('\n').rstrip(';')

            # Remove leading and trailing double quotes.
            if ligne.startswith('"'):
                ligne = ligne[1:]
            if ligne.endswith('"'):
                ligne = ligne[:-1]

            # Replace double quotes with single quotes inside cells.
            ligne = ligne.replace('""', '"')

            lignes_propres.append(ligne)

    # Keep only lines with the expected number of columns.
    # Malformed lines (108 total) correspond to Probe entries
    # that would be removed in step 1 anyway.
    lignes_valides = []
    for ligne in lignes_propres:
        nb = len(list(csv.reader([ligne]))[0])
        if nb == EXPECTED_COLUMNS:
            lignes_valides.append(ligne)

    # Reassemble lines and parse with csv.reader.
    contenu_propre = '\n'.join(lignes_valides)
    reader = csv.reader(StringIO(contenu_propre))
    lignes = [ligne for ligne in reader]

    colonnes = lignes[0]
    donnees = lignes[1:]

    df = pd.DataFrame(donnees, columns=colonnes)

    # Strip invisible whitespace from all text columns.
    df = df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)

    print(f"{nom_fichier} chargé : {df.shape[0]} lignes, {df.shape[1]} colonnes")
    return df


# ============================================================
# ÉTAPE 1 : SUPPRESSION DES DONNÉES AMBIGUËS
# ============================================================

def etape1(df, nom):
    """Remove ambiguous observations from the dataset.

    Two types of rows are removed:
        - Rows where Activity == 'Probe' (biological probes,
          not real inhibitor candidates).
        - Rows where Activity == 'Active' but Activity_Value is
          empty (cannot be verified against the 10 µM threshold).

    Parameters
    ----------
    df : pandas.DataFrame
        Raw input DataFrame.
    nom : str
        Protein name, used for display purposes.

    Returns
    -------
    pandas.DataFrame
        Cleaned DataFrame without ambiguous observations.
    """
    n0 = len(df)

    # Remove Probe entries.
    df = df[df["Activity"] != "Probe"]
    n1 = len(df)

    # Remove Active entries without a numerical value.
    masque = (df["Activity"] == "Active") & (df["Activity_Value"] == "")
    df = df[~masque]
    n2 = len(df)

    print(
        f"{nom} | Avant : {n0} | "
        f"Après suppression Probe : {n1} | "
        f"Après suppression Active sans valeur : {n2}"
    )
    return df


# ============================================================
# ÉTAPE 2a : RÉÉTIQUETAGE SELON LE SEUIL 10 µM
# ============================================================

def etape2a(df, nom):
    """Relabel activity based on the 10 µM threshold.

    Each observation with a numerical Activity_Value is reclassified:
        - Activity_Value <= 10 µM  ->  'Active'
        - Activity_Value >  10 µM  ->  'Inactive'
        - No Activity_Value        ->  label kept as is

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame after step 1.
    nom : str
        Protein name, used for display purposes.

    Returns
    -------
    pandas.DataFrame
        DataFrame with updated Activity column.
    """
    # Convert Activity_Value to float; non-convertible values become NaN.
    df["Activity_Value"] = pd.to_numeric(df["Activity_Value"], errors="coerce")

    def requalify(row):
        """Relabel a single row based on the activity threshold."""
        if pd.notna(row["Activity_Value"]):
            if row["Activity_Value"] <= ACTIVITY_THRESHOLD:
                return "Active"
            else:
                return "Inactive"
        return row["Activity"]

    df["Activity"] = df.apply(requalify, axis=1)

    counts = df["Activity"].value_counts()
    print(f"{nom} | {counts.to_dict()}")
    return df


# ============================================================
# ÉTAPE 2b : CONSOLIDATION PAR MOLÉCULE (SID)
# ============================================================

def consolider_groupe(groupe):
    """Apply the decision tree to consolidate a group of rows for one SID.

    A molecule (identified by its SID) can appear multiple times
    if tested in several experiments. This function merges all
    observations into a single row using the following rules:

        - Case 1: All labels identical -> keep label, average values.
        - Case 2a: Contradictory labels, all values present:
            * mean > std  -> keep molecule, label by threshold.
            * mean <= std -> eliminate (too variable).
        - Case 2b: Contradictory labels, some values missing:
            * Clear majority label -> keep majority label.
            * No majority          -> eliminate.

    Parameters
    ----------
    groupe : pandas.DataFrame
        Group of rows sharing the same Substance_SID.

    Returns
    -------
    pandas.Series
        A Series with Activity, Activity_Value and a boolean 'garder'
        indicating whether the molecule should be kept.
    """
    labels = groupe["Activity"].tolist()
    valeurs = groupe["Activity_Value"].dropna().tolist()
    unique_labels = set(labels)

    # Case 1: all labels are identical.
    if len(unique_labels) == 1:
        label = labels[0]
        valeur = sum(valeurs) / len(valeurs) if valeurs else None
        return pd.Series({"Activity": label, "Activity_Value": valeur, "garder": True})

    # Case 2: contradictory labels.
    else:
        toutes_valeurs = groupe["Activity_Value"].isna().sum() == 0

        # Sub-case 2a: all rows have a numerical value.
        if toutes_valeurs:
            moyenne = sum(valeurs) / len(valeurs)
            ecart_type = pd.Series(valeurs).std()
            if moyenne > ecart_type:
                label = "Active" if moyenne <= ACTIVITY_THRESHOLD else "Inactive"
                return pd.Series({"Activity": label, "Activity_Value": moyenne, "garder": True})
            else:
                # Too much variability: eliminate the molecule.
                return pd.Series({"Activity": None, "Activity_Value": None, "garder": False})

        # Sub-case 2b: at least one row has no numerical value.
        else:
            compte = Counter(labels)
            max_compte = max(compte.values())
            majorite = [k for k, v in compte.items() if v == max_compte]

            if len(majorite) == 1:
                label = majorite[0]
                valeur = sum(valeurs) / len(valeurs) if (label == "Active" and valeurs) else None
                return pd.Series({"Activity": label, "Activity_Value": valeur, "garder": True})
            else:
                # No clear majority: eliminate the molecule.
                return pd.Series({"Activity": None, "Activity_Value": None, "garder": False})


def etape2b(df, nom):
    """Consolidate multiple observations per molecule into one row.

    Groups rows by Substance_SID and applies the decision tree
    defined in consolider_groupe(). Metadata (CID, Activity_Type,
    Activity_Qualifier) is preserved from the first occurrence.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame after step 2a.
    nom : str
        Protein name, used for display purposes.

    Returns
    -------
    pandas.DataFrame
        DataFrame with one row per unique molecule (SID).
    """
    # Keep metadata from the first occurrence of each SID.
    meta = df.drop_duplicates(subset=["Substance_SID"])[
        ["Substance_SID", "Compound_CID", "Activity_Type", "Activity_Qualifier"]
    ].copy()

    # Apply the decision tree to each group of rows sharing the same SID.
    consolide = df.groupby("Substance_SID").apply(
        consolider_groupe, include_groups=False
    ).reset_index()

    # Remove molecules flagged for elimination.
    consolide = consolide[consolide["garder"] == True].drop(columns=["garder"])

    # Merge metadata with consolidated activity data.
    result = meta.merge(consolide, on="Substance_SID", how="inner")

    print(f"{nom} | Molécules uniques : {len(result)} | {result['Activity'].value_counts().to_dict()}")
    return result


# ============================================================
# ÉTAPE 3 : SUPPRESSION DES FAUX POSITIFS
# ============================================================

def get_cids_actifs(aid):
    """Retrieve the set of active CIDs for a given PubChem assay (AID).

    Parameters
    ----------
    aid : int
        PubChem Assay ID (AID).

    Returns
    -------
    set of str
        Set of active Compound CIDs as strings.
        Returns an empty set if the request fails.
    """
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        f"/assay/aid/{aid}/cids/JSON?active=true"
    )
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        data = r.json()
        cids = data["InformationList"]["Information"][0].get("CID", [])
        cids = set(str(c) for c in cids)
        print(f"  AID {aid} : {len(cids)} CIDs actifs récupérés")
        return cids
    except Exception as e:
        print(f"  AID {aid} : erreur ({e})")
        return set()


def get_faux_positifs():
    """Retrieve all CIDs identified as potential false positives.

    Three categories of false positives are handled:
        1. Luciferase inhibitors: active in counter-screens (AID 585,
           485341) but absent from confirm-screens (AID 584, 485294).
        2. Promiscuous aggregators: active in AID 411.
        3. Auto-fluorescent compounds: active in AID 587-594.

    Returns
    -------
    set of str
        Union of all false positive CIDs as strings.
    """
    print("Récupération des faux positifs depuis PubChem...")

    # Category 1: luciferase inhibitors.
    cids_counter = get_cids_actifs(585) | get_cids_actifs(485341)
    cids_confirm = get_cids_actifs(584) | get_cids_actifs(485294)
    fp_luciferase = cids_counter - cids_confirm
    print(f"  Faux positifs luciférase : {len(fp_luciferase)}")

    # Category 2: promiscuous aggregators.
    fp_agregateurs = get_cids_actifs(411)
    print(f"  Faux positifs agrégateurs : {len(fp_agregateurs)}")

    # Category 3: auto-fluorescent compounds.
    fp_fluor = set()
    for aid in AIDS_AUTOFLUORESCENT:
        fp_fluor |= get_cids_actifs(aid)
    print(f"  Faux positifs auto-fluorescents : {len(fp_fluor)}")

    tous = fp_luciferase | fp_agregateurs | fp_fluor
    print(f"  Total faux positifs CIDs : {len(tous)}")
    return tous


def etape3(df, fp_cids, nom):
    """Remove active molecules identified as potential false positives.

    Active molecules whose Compound_CID is found in the false positive
    list are removed from the dataset.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame after step 2b.
    fp_cids : set of str
        Set of Compound CIDs identified as false positives.
    nom : str
        Protein name, used for display purposes.

    Returns
    -------
    pandas.DataFrame
        DataFrame without false positive active molecules.
    """
    n0 = len(df)

    # Strip whitespace from CIDs before comparison.
    df["Compound_CID"] = df["Compound_CID"].str.strip()

    masque = (df["Activity"] == "Active") & (df["Compound_CID"].isin(fp_cids))
    df = df[~masque]

    n1 = len(df)
    print(f"{nom} | Faux positifs supprimés : {n0 - n1} | Restant : {n1}")
    return df


# ============================================================
# ÉTAPE 4 : EXPORT DES CSV FINAUX
# ============================================================

def etape4(df, nom):
    """Export the final cleaned dataset to a CSV file.

    Only the 6 columns required by the project are kept,
    in the exact order specified in the assignment.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame after step 3.
    nom : str
        Protein name, used to build the output filename.

    Returns
    -------
    pandas.DataFrame
        Final DataFrame with only the 6 required columns.
    """
    df_final = df[FINAL_COLUMNS].copy()

    nom_fichier = f"{nom}_final.csv"
    df_final.to_csv(nom_fichier, index=False)

    print(f"{nom} | Fichier exporté : {nom_fichier} | {len(df_final)} lignes")
    print(f"  Actifs : {(df_final['Activity'] == 'Active').sum()}")
    print(f"  Inactifs : {(df_final['Activity'] == 'Inactive').sum()}")
    print(f"  Unspecified : {(df_final['Activity'] == 'Unspecified').sum()}")
    return df_final


# ============================================================
# ÉTAPE 5 : ANALYSE CROISÉE
# ============================================================

def etape5(df_b1, df_g2, df_c1):
    """Compute cross-protein analysis of active molecules.

    Identifies molecules active for exactly one, two or all three
    proteins, as well as molecules never active for any protein,
    using set operations on Compound CIDs.

    Parameters
    ----------
    df_b1 : pandas.DataFrame
        Final DataFrame for ABCB1.
    df_g2 : pandas.DataFrame
        Final DataFrame for ABCG2.
    df_c1 : pandas.DataFrame
        Final DataFrame for ABCC1.

    Returns
    -------
    dict
        Dictionary containing sets for each category:
        seulement_b1, seulement_g2, seulement_c1,
        deux_b1_g2, deux_b1_c1, deux_g2_c1,
        trois, jamais_actifs.
    """
    actifs_b1 = set(df_b1[df_b1["Activity"] == "Active"]["Compound_CID"])
    actifs_g2 = set(df_g2[df_g2["Activity"] == "Active"]["Compound_CID"])
    actifs_c1 = set(df_c1[df_c1["Activity"] == "Active"]["Compound_CID"])

    tous_cids = (
        set(df_b1["Compound_CID"])
        | set(df_g2["Compound_CID"])
        | set(df_c1["Compound_CID"])
    )

    # Molecules active for exactly one protein.
    seulement_b1 = actifs_b1 - actifs_g2 - actifs_c1
    seulement_g2 = actifs_g2 - actifs_b1 - actifs_c1
    seulement_c1 = actifs_c1 - actifs_b1 - actifs_g2

    # Molecules active for exactly two proteins.
    deux_b1_g2 = (actifs_b1 & actifs_g2) - actifs_c1
    deux_b1_c1 = (actifs_b1 & actifs_c1) - actifs_g2
    deux_g2_c1 = (actifs_g2 & actifs_c1) - actifs_b1

    # Pan-inhibitors: active for all three proteins.
    trois = actifs_b1 & actifs_g2 & actifs_c1

    # Molecules never active for any protein.
    jamais_actifs = tous_cids - actifs_b1 - actifs_g2 - actifs_c1

    print(f"Actifs seulement ABCB1          : {len(seulement_b1)}")
    print(f"Actifs seulement ABCG2          : {len(seulement_g2)}")
    print(f"Actifs seulement ABCC1          : {len(seulement_c1)}")
    print(f"Actifs ABCB1 + ABCG2            : {len(deux_b1_g2)}")
    print(f"Actifs ABCB1 + ABCC1            : {len(deux_b1_c1)}")
    print(f"Actifs ABCG2 + ABCC1            : {len(deux_g2_c1)}")
    print(f"Actifs pour les 3 protéines     : {len(trois)}")
    print(f"Jamais actifs pour aucune       : {len(jamais_actifs)}")

    return {
        "seulement_b1": seulement_b1,
        "seulement_g2": seulement_g2,
        "seulement_c1": seulement_c1,
        "deux_b1_g2": deux_b1_g2,
        "deux_b1_c1": deux_b1_c1,
        "deux_g2_c1": deux_g2_c1,
        "trois": trois,
        "jamais_actifs": jamais_actifs,
    }


# ============================================================
# ÉTAPE 6 : VISUALISATIONS
# ============================================================

def etape6(dfs_finaux, croise):
    """Generate and save visualisation figures.

    Produces two PNG figures:
        1. Bar chart of activity distribution per protein
           (Unspecified excluded for readability).
        2. Bar chart of active molecule overlaps between proteins
           (Never active excluded for readability).

    Parameters
    ----------
    dfs_finaux : list of pandas.DataFrame
        List of final DataFrames [ABCB1, ABCG2, ABCC1].
    croise : dict
        Dictionary returned by etape5() with overlap sets.

    Returns
    -------
    None
    """
    proteines = ["ABCB1", "ABCG2", "ABCC1"]
    couleurs = {
        "Active": "#2ecc71",
        "Inactive": "#e74c3c",
        "Unspecified": "#95a5a6",
        "Inconclusive": "#f39c12",
    }

    # Figure 1: activity distribution without Unspecified.
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.suptitle(
        "Distribution des molécules par activité\n(sans Unspecified)",
        fontsize=14,
        fontweight="bold",
    )

    for ax, nom, df in zip(axes, proteines, dfs_finaux):
        df_filtre = df[df["Activity"] != "Unspecified"]
        counts = df_filtre["Activity"].value_counts()
        ax.bar(
            counts.index,
            counts.values,
            color=[couleurs.get(l, "gray") for l in counts.index],
            edgecolor="white",
            width=0.5,
        )
        ax.set_title(nom, fontsize=12, fontweight="bold")
        ax.set_ylabel("Nombre de molécules")
        for i, val in enumerate(counts.values):
            ax.text(i, val + 2, str(val), ha="center", fontweight="bold")

        n_unspec = (df["Activity"] == "Unspecified").sum()
        ax.set_xlabel(f"(+ {n_unspec} Unspecified non affichés)", fontsize=8)

    plt.tight_layout()
    plt.savefig("distribution_activite.png", dpi=150)
    print("Figure sauvegardée : distribution_activite.png")
    plt.close()

    # Figure 2: active molecule overlaps without Never active.
    fig, ax = plt.subplots(figsize=(12, 6))
    categories = [
        "Seul.\nABCB1", "Seul.\nABCG2", "Seul.\nABCC1",
        "ABCB1+\nABCG2", "ABCB1+\nABCC1", "ABCG2+\nABCC1",
        "Les 3\nprotéines",
    ]
    valeurs = [
        len(croise["seulement_b1"]),
        len(croise["seulement_g2"]),
        len(croise["seulement_c1"]),
        len(croise["deux_b1_g2"]),
        len(croise["deux_b1_c1"]),
        len(croise["deux_g2_c1"]),
        len(croise["trois"]),
    ]
    palette = [
        "#3498db", "#9b59b6", "#e67e22",
        "#1abc9c", "#f39c12", "#e74c3c", "#2c3e50",
    ]

    bars = ax.bar(categories, valeurs, color=palette, edgecolor="white", linewidth=1.5)
    ax.set_title(
        f"Chevauchements entre molécules actives\n"
        f"(+ {len(croise['jamais_actifs'])} jamais actifs non affichés)",
        fontsize=13,
        fontweight="bold",
    )
    ax.set_ylabel("Nombre de molécules")
    for bar, val in zip(bars, valeurs):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 2,
            str(val),
            ha="center",
            fontweight="bold",
        )

    plt.tight_layout()
    plt.savefig("chevauchements.png", dpi=150)
    print("Figure sauvegardée : chevauchements.png")
    plt.close()


# ============================================================
# PROGRAMME PRINCIPAL
# ============================================================

if __name__ == "__main__":

    print("=== CHARGEMENT DES FICHIERS ===")
    df_abcb1 = charger_csv("ABCB1.csv")
    df_abcg2 = charger_csv("ABCG2.csv")
    df_abcc1 = charger_csv("ABCC1.csv")

    print("\n=== ÉTAPE 1 : SUPPRESSION DES DONNÉES AMBIGUËS ===")
    df_abcb1 = etape1(df_abcb1, "ABCB1")
    df_abcg2 = etape1(df_abcg2, "ABCG2")
    df_abcc1 = etape1(df_abcc1, "ABCC1")

    print("\n=== ÉTAPE 2a : RÉÉTIQUETAGE ===")
    df_abcb1 = etape2a(df_abcb1, "ABCB1")
    df_abcg2 = etape2a(df_abcg2, "ABCG2")
    df_abcc1 = etape2a(df_abcc1, "ABCC1")

    print("\n=== ÉTAPE 2b : CONSOLIDATION PAR MOLÉCULE ===")
    df_abcb1 = etape2b(df_abcb1, "ABCB1")
    df_abcg2 = etape2b(df_abcg2, "ABCG2")
    df_abcc1 = etape2b(df_abcc1, "ABCC1")

    print("\n=== ÉTAPE 3 : SUPPRESSION DES FAUX POSITIFS ===")
    fp_cids = get_faux_positifs()
    df_abcb1 = etape3(df_abcb1, fp_cids, "ABCB1")
    df_abcg2 = etape3(df_abcg2, fp_cids, "ABCG2")
    df_abcc1 = etape3(df_abcc1, fp_cids, "ABCC1")

    print("\n=== ÉTAPE 4 : EXPORT DES CSV FINAUX ===")
    df_abcb1_final = etape4(df_abcb1, "ABCB1")
    df_abcg2_final = etape4(df_abcg2, "ABCG2")
    df_abcc1_final = etape4(df_abcc1, "ABCC1")

    print("\n=== ÉTAPE 5 : ANALYSE CROISÉE ===")
    croise = etape5(df_abcb1_final, df_abcg2_final, df_abcc1_final)

    print("\n=== ÉTAPE 6 : VISUALISATIONS ===")
    etape6([df_abcb1_final, df_abcg2_final, df_abcc1_final], croise)

    print("\n=== ANALYSE TERMINÉE ===")
