import pandas as pd
from datetime import date

mut_df = pd.read_csv("results/mutation_frequencies_per_position.csv")
var_df = pd.read_csv("results/variant_frequencies.csv")

n_seqs = 250
n_mut_positions = len(mut_df)
n_high_freq = len(mut_df[mut_df["mutation_frequency"] >= 0.05])
most_mutated = mut_df.loc[mut_df["mutation_frequency"].idxmax()]

report = f"""# blaNDM Analysis Summary
**Date**: {date.today().isoformat()}
**Dataset**: {n_seqs} blaNDM sequences from *Klebsiella pneumoniae*

## Key Findings
- **Mutation analysis:**
  - Total variable positions: {n_mut_positions}
  - Highly variable positions (>=5% frequency): {n_high_freq}
  - Most mutated position: pos {int(most_mutated['position'])} ({most_mutated['mutation_frequency']*100:.1f}% of sequences)

## Variant Distribution
Top 10 most common patterns:
{var_df.head(10).to_string(index=False)}

Full table: `results/variant_frequencies.csv`

## Phylogeny
- **Input**: `data/blaNDM_Kpneumoniae_phylogeny.fasta` (251 sequences)
- **Tree file**: `data/blaNDM_Kpneumoniae_phylogeny.fasta.contree`
- **Visualization**: `figures/phylogeny_tree.svg` and `figures/phylogeny_tree.html`
- **Combined figure**: `figures/combined_figure.png`

## Output Files
- `figures/combined_figure.png` -- Tree + mutation heatmap
- `figures/phylogeny_tree.svg` -- Tree visualization (vector)
- `results/mutation_frequencies_per_position.csv` -- Per-position mutation data
- `results/variant_frequencies.csv` -- Variant pattern frequencies
"""

with open("results/analysis_summary.md", "w") as f:
    f.write(report)
print(report)
