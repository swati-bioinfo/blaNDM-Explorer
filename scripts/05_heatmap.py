import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from pathlib import Path

DATA = Path("data")
RESULTS = Path("data")

aln = list(SeqIO.parse(DATA / "blaNDM_Kpneumoniae_aligned_trimmed.fasta", "fasta"))

try:
    meta = pd.read_csv(DATA / "sequence_metadata.csv")
    meta["year_num"] = pd.to_numeric(meta["year"], errors="coerce")
    meta = meta.sort_values("year_num", na_position="last")
    ordered_ids = meta["short_id"].tolist()
except:
    ordered_ids = [rec.id for rec in aln]

df_muts = pd.read_csv(RESULTS / "all_mutations.csv")

pivot = df_muts.pivot_table(index="sequence", columns="position",
                           values="is_mutation", aggfunc="max", fill_value=False)

pivot = pivot.reindex([i for i in ordered_ids if i in pivot.index]).astype(float)

variable_cols = [c for c in pivot.columns if pivot[c].sum() > 0]
pivot = pivot[variable_cols]

fig, ax = plt.subplots(figsize=(14, max(8, len(pivot) * 0.15)))
sns.heatmap(pivot, cmap="YlOrRd", cbar_kws={"label": "Mutation present"},
            linewidths=0, ax=ax, vmin=0, vmax=1)
ax.set_xlabel("Position (bp)", fontsize=12)
ax.set_ylabel("Sequences (sorted by year)", fontsize=12)
ax.set_title("Mutation Heatmap of blaNDM — K. pneumoniae", fontsize=14)

plt.tight_layout()
plt.savefig(RESULTS / "mutation_heatmap.png", dpi=300)
plt.show()
print(f"Saved: {RESULTS / 'mutation_heatmap.png'} (variable positions: {len(variable_cols)})")
