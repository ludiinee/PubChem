# ABC Efflux Pump Inhibitor Analysis

Analyse bioinformatique des données expérimentales issues de PubChem pour trois transporteurs ABC efflux (ABCB1, ABCG2, ABCC1). Le script identifie leurs inhibiteurs potentiels à travers un pipeline de nettoyage et d'analyse multi-étapes.

**Auteure :** Ludine  
**Encadrant :** Dr. TRAN NGUYEN Viet Khoa  
**Année :** 2026  
**Version :** 1.0.0

---

## Contexte

Les transporteurs ABC (ATP-Binding Cassette) jouent un rôle central dans la résistance aux médicaments. Ce projet analyse des données de criblage haut-débit (HTS) issues de PubChem afin d'identifier des molécules inhibant les protéines ABCB1, ABCG2 et ABCC1, avec un seuil d'activité fixé à **10 µM**.

---

## Prérequis

- Python 3.8+
- Les bibliothèques suivantes :
```bash
pip install pandas matplotlib requests
```

---

## Fichiers d'entrée

Les trois fichiers CSV suivants doivent être présents dans le **même dossier** que le script :

| Fichier | Protéine ciblée |
|---|---|
| `ABCB1.csv` | P-glycoprotéine (MDR1) |
| `ABCG2.csv` | Breast Cancer Resistance Protein |
| `ABCC1.csv` | Multidrug Resistance Protein 1 |

Ces fichiers sont téléchargeables depuis [PubChem BioAssay](https://pubchem.ncbi.nlm.nih.gov/bioassay/).

---

## Utilisation
```bash
python analyse2.py
```

---

## Pipeline d'analyse

| Étape | Description |
|---|---|
| **Étape 0** | Chargement et nettoyage du format CSV PubChem |
| **Étape 1** | Suppression des entrées ambiguës (Probe, Active sans valeur) |
| **Étape 2a** | Réétiquetage selon le seuil 10 µM |
| **Étape 2b** | Consolidation par molécule (SID) via arbre de décision |
| **Étape 3** | Suppression des faux positifs via l'API PubChem |
| **Étape 4** | Export des CSV finaux |
| **Étape 5** | Analyse croisée inter-protéines |
| **Étape 6** | Visualisations (graphiques PNG) |

---

## Fichiers de sortie

| Fichier | Contenu |
|---|---|
| `ABCB1_final.csv` | Données nettoyées pour ABCB1 |
| `ABCG2_final.csv` | Données nettoyées pour ABCG2 |
| `ABCC1_final.csv` | Données nettoyées pour ABCC1 |
| `distribution_activite.png` | Distribution Active/Inactive par protéine |
| `chevauchements.png` | Chevauchements des molécules actives entre protéines |

Chaque CSV final contient 6 colonnes : `Substance_SID`, `Compound_CID`, `Activity`, `Activity_Type`, `Activity_Qualifier`, `Activity_Value`.

---

## Notes techniques

- Le format des CSV PubChem nécessite un prétraitement ligne par ligne avant parsing standard (guillemets imbriqués, points-virgules en fin de ligne, champs multi-lignes).
- La détection des faux positifs utilise le **Compound_CID** plutôt que le Substance_SID, car les SIDs proviennent de ChEMBL et ne correspondent pas aux SIDs natifs retournés par l'API PubChem.
- Les AIDs utilisés pour la détection des faux positifs : luciférase (585, 485341, 584, 485294), agrégateurs (411), autofluorescence (587–594).
