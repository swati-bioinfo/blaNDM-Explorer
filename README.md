# blaNDM Explorer
Interactive phylogenetic and mutation analysis of 250 NDM carbapenemase sequences from *Klebsiella pneumoniae*.
**Live app**: [https://blandm-explorer.streamlit.app](https://blandm-explorer.streamlit.app)
---
## Overview
This project performs a complete bioinformatics pipeline on NDM metallo-beta-lactamase variants:
- **250 sequences** downloaded from NCBI Virus (blaNDM gene, *K. pneumoniae*)
- **Multiple sequence alignment** via MAFFT, trimmed with ClipKIT
- **Maximum Likelihood phylogeny** with IQ-TREE (HKY+G4, 1000 bootstrap)
- **Mutation frequency analysis** against NDM-1 reference (FN396876)
- **Interactive visualization** via Streamlit + toytree
## Key Findings
| Metric | Value |
|--------|-------|
| Alignment length | 2,271 bp |
| Parsimony-informative sites | 1,343 |
| Unique sequence types | 90 (after 161 duplicates removed) |
| NDM variant patterns | 63 |
| Best-fit model | HKY+G4 |
| Log-likelihood | -15,943.42 |
| Bootstrap correlation | 0.990 |
| Insertion hotspot | Positions 789-795 (90.4% of sequences) |
| Most divergent taxon | NDM-63 (~0.68 subs/site from others) |
## Repository Structure
├── scripts/
│   ├── 01_rename_sequences.py      # Rename NCBI headers → short IDs
│   ├── 02_alignment_report.py      # Alignment statistics
│   ├── 02b_combine_reference.py    # Add NDM-1 reference
│   ├── 03_mutation_frequency.py    # Per-position mutation analysis
│   ├── 04_barplot.py               # Mutation frequency bar chart
│   ├── 05_heatmap.py               # Mutation heatmap
│   ├── 06_variant_table.py         # Variant frequency table
│   ├── 07_deduplicate.py           # Remove identical sequences
│   ├── 08_rename_for_phylogeny.py  # NDM-variant_Accession labels
│   ├── 09_tree_visualize.py        # Toytree tree visualization
│   ├── 10_combine_figure.py        # Tree + heatmap composite
│   ├── 11_summary_report.py        # Final report generation
│   ├── streamlit_app.py            # Interactive Streamlit app
│   └── requirements.txt            # Python dependencies
├── figures/                        # Generated visualizations
├── data/                           # Sequence files (FASTA, CSV)
├── results/                        # Analysis outputs (CSV, PNG)
└── .gitignore
## Interactive App Pages
| Page | Description |
|------|-------------|
| **Overview** | 8 KPI metrics, IQ-TREE summary, key findings, saved plots |
| **Interactive Tree** | Full phylogeny rendered with toytree, colored by NDM variant (29 colors), exportable HTML |
| **Mutation Analysis** | Per-position frequency with threshold slider, position explorer, saved heatmap + barplot |
| **Variant Distribution** | Sorted bar chart + table of all 63 variant patterns |
| **Combined Figure** | Tree + heatmap composite side-by-side |
## Pipeline
1. **Download** 250 blaNDM sequences from NCBI Virus
2. **Rename** headers to clean IDs (`KpNDM_0001`)
3. **Align** with MAFFT, trim with ClipKIT
4. **Add** NDM-1 reference (FN396876), re-align
5. **Compute** per-position mutation frequencies
6. **Rename** for phylogeny (`NDM-63_OR667647.1`)
7. **Build** ML tree with IQ-TREE (HKY+G4, -bb 1000)
8. **Visualize** with toytree (tip colors by variant, bootstrap labels)
## Local Setup
```bash
# Clone the repo
git clone https://github.com/swati-bioinfo/blaNDM-Explorer.git
cd blaNDM-Explorer
# Install dependencies
pip install -r scripts/requirements.txt
# Run the analysis (each step)
py scripts/01_rename_sequences.py
py scripts/03_mutation_frequency.py
# ... etc.
# Launch the interactive app
py -m streamlit run scripts/streamlit_app.py
Deployment
The app is deployed on Streamlit Cloud (https://share.streamlit.io). Push to GitHub triggers automatic redeployment:
git add .
git commit -m "Update"
git push
Dependencies
- Python 3.11+
- streamlit, pandas, toytree, toyplot
- Biopython (pipeline scripts)
- ClipKIT (alignment trimming)
License
Educational / research use.
