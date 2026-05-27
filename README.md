# splice-induced-ppi-loss
Computational analysis of splice-induced disruption of protein-protein interactions and their structural interfaces
# Splice-Induced Loss of Protein-Protein Interactions and Impact on Protein Interaction Interfaces

Bachelor's thesis — Vrije Universiteit Amsterdam  
Computational analysis of how alternative splicing disrupts protein-protein interactions , using PiPENN interface predictions, pairwise sequence alignment, and statistical testing across 8 deep learning models.

---

## Project Overview

Alternative splicing can remove regions of a protein that are critical for binding to interaction partners. This project investigates whether residues spliced out in alternative isoforms are disproportionately located at predicted protein-protein interaction (PPI) interfaces, and whether this enrichment correlates with experimentally observed loss of interaction.

The pipeline:
1. Filters isoforms and interactions from the Fixed Interactor Panel
2. Aligns alternative isoforms to reference sequences to identify spliced-out positions
3. Overlays spliced positions with PiPENN interface predictions
4. Tests whether spliced residues are enriched at predicted interfaces
5. Maps spliced interface residues onto Pfam domain annotations
6. Applies Benjamini-Hochberg correction across 8 models

---

## Repository Structure

```
.
│
│── orfs_translation.py                                # Translate ORF nucleotide sequences to protein (stop at first stop codon)
├── Filtration.py                                      # Filter interactions: remove N/A rows, keep canonical-positive pairs only
├── Biopyhton_pairwise_alignment.py                    # Global pairwise alignment of isoforms vs. reference; identify spliced-out positions
├── Final_table.py                                     # Integrate PiPENN scores with spliced positions into master residue table
│
│   # -- Statistical analysis --
├── MO-GI.py                                           # Identify mixed-outcome gene-interactor pairs (one positive + one negative isoform)
├── Statistical_analysis.py                            # Wilcoxon, Jaccard, and chi-square tests per model; append p-values to all_pvalues.csv
├── BH.py                                              # Benjamini-Hochberg FDR correction across all 8 models
│
│   # -- Exploratory analysis --
├── Step_1_Unique_Canonical_Sequences_Filtration.py    # Extract canonical sequences for genes present in the dataset
├──Step_2_Interproscan_Rest_API.p                      # fetches Pfam domain annotations via InterPro REST API;submits canonical sequences
├── Step_3_and_4.py                                    # Expand to residue level; annotate with Pfam domains; calculate prop_A/B/C
├── Step_5_domain_plot.py                              # Generates lost vs. preserved side-by-side domain overlap bar charts (PDF)
├── One_Histogram_formation.py                         # Histogram: percentage of predicted interface lost per isoform
└── Density_plot.py                                    # Overlaid KDE density plots of interface residues lost across all 8 models
```

---

## Pipeline

Two input streams are processed in parallel before merging:

```
SPLICED PROTEIN SIDE                    INTERACTOR SIDE
────────────────────                    ───────────────
Isoform ORF dataset                     Fixed interactor panel
(1,423 isoforms · Yang et al. 2016)     (502 proteins · ORFeome 7.1)
        ↓                                       ↓
orfs_translation.py                     Filtration.py
(ORF → protein, stop at first codon)    (remove N/A rows; keep
        ↓                                canonical-positive pairs only)
PiPENN-EMB interface prediction                 ↓
(6 local models + Web-DNET)                     │
per-residue interface scores                    │
        ↓                                       │
        └──────────────┬─────────────────────────┘
                       ↓
           Biopyhton_pairwise_alignment.py
           (isoform vs. canonical protein · global mode)
           (primary: isoform _1 as reference)
           (fallback: longest isoform · 0 spliced positions excluded)
                       ↓
               Final_table.py
               (data integration: master residue table)
               (join key: Gene_Symbol / Isoform_ID)
                       ↓
        ┌──────────────┼──────────────────────┐
        ↓              ↓                      ↓
Jaccard index    Chi-squared test     Wilcoxon rank-sum (PRIMARY)
(interface ∩     (2×2 residue-level   (unpaired · lost vs. preserved
 spliced /        contingency table)   interaction group)
 interface ∪                                  ↓
 spliced)                            MO-GI.py
        │              │             (mixed-outcome gene–interactor pairs)
        │              │                      ↓
        │              │             Statistical_analysis.py
        │              │             (Jaccard · χ² · Wilcoxon repeated
        │              │              on MO-GI dataset, same parameters)
        └──────────────┴──────────────────────┘
                       ↓
                    BH.py
          (Benjamini-Hochberg FDR correction
           across all 8 models · all_pvalues.csv)
                       ↓
                 Interpretation
    Does splicing remove high-scoring interface residues,
         leading to loss of protein–protein interactions?
```

**Exploratory analysis** (run after `Final_table.py`, no correction applied):
```
Step_1_Unique_Canonical_Sequences_Filtration.py → extract canonical sequences
Step_2_Interproscan_Rest_API.py                 → submits canonical sequences
Step_3_and_4.py                                 → residue expansion + Pfam annotation + prop_A/B/C
Step_5_plot_formation.py                        → domain overlap bar charts lost vs preserved (PDF, 24 plots)
One_Histogram_formation.py                      → histogram: % of interface lost per isoform
Density_plot.py                                 → overlaid KDE density plots across all 8 models
```

---

## Input Files Required

| File | Description |
|------|-------------|
| `preds-ensemble-all.csv` | PiPENN predictions (ensemble model) |
| `preds-<model>-all.csv` | PiPENN predictions for each of the 8 models |
| `Fixed_interactor_table.xlsx` | Interactor panel with interaction outcomes (Yang et al.) |
| `protein_sequences.fasta` | All isoform protein sequences |
| `ORFs_sequences.fasta` | ORF nucleotide sequences (for `orfs_translation.py`) |


---

## Output Files

| File | Description |
|------|-------------|
| `Filtered_interactor_table.csv` | Interactions passing all 3 filtration criteria |
| `spliced_positions_reference1.csv` | Spliced positions using isoform _1 as reference |
| `spliced_positions_longest.csv` | Spliced positions using longest isoform as reference |
| `preds-<model>-all_final_table_reference1.csv` | Final table with PiPENN scores per isoform |
| `Triples_<Model>.csv` | Mixed-outcome pairs (one positive + one negative isoform per gene-interactor pair) |
| `interpro_pfam_only.tsv` | Pfam domain annotations (from InterPro) |
| `domain_analysis_results/` | Per-model proportion CSVs (prop_A, prop_B, prop_C) |
| `all_pvalues.csv` | Raw p-values across all 8 models |
| `BH_corrected_results.xlsx` | BH-corrected p-values and significance markers |
| `all_24_plots_lost_vs_preserved.pdf` | Domain overlap bar charts lost vs preserved for all models and groups |
| `interface_lost_all_models_overlaid.png` | Overlaid KDE density plot |

---

## Methods Summary

**Interface prediction** — PiPENN (residue-level scores; threshold > 0.5 = interface)

**Alignment** — Global pairwise alignment via BioPython `PairwiseAligner` (match +1, mismatch −1, gap open −10, gap extend −0.5); spliced positions = reference residues absent from alternative isoform

**Statistical tests (per model):**
- Wilcoxon rank-sum test — mean PiPENN scores at spliced positions, negative vs positive interaction outcome
- Jaccard index + Wilcoxon — overlap between spliced positions and predicted interface, negative vs positive
- Chi-square test — residue-level contingency table (spliced × interface)

**Multiple testing correction** — Benjamini-Hochberg FDR across 8 models (`BH.py`)

**Domain analysis** — Pfam domain annotations from InterPro; proportions calculated per isoform × domain:
- `prop_A` — fraction of domain predicted as interface
- `prop_B` — fraction of domain that is spliced out
- `prop_C` — fraction of domain that is both spliced out and predicted as interface (key metric)

---

## Models

Eight PiPENN model variants were evaluated in parallel:

| Label | Model |
|-------|-------|
| Ensemble | ensemble |
| RNET | rnet |
| CNN-RNN | cnn-rnn |
| DNET | dnet |
| ANN | ann |
| UNET | unet |
| RNN | rnn |
| Web-DNET | web-dnet |

---

## Dependencies

```bash
pip install biopython pandas scipy numpy matplotlib openpyxl
```

Python 3.8+

---

## Notes

- All scripts use the `reference1` (isoform `_1`) strategy by default; `longest` reference variants are also generated in `Biopyhton_pairwise_alignment.py`
- `Statistical_analysis.py` is designed to be run once per model file; it appends one row to `all_pvalues.csv` each run, which `BH.py` then reads
- Gain-of-sequence isoforms (no spliced-out positions after alignment) are excluded from all downstream analyses
