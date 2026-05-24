import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

DATA = Path("data")

df = pd.read_csv(DATA / "mutation_frequencies_per_position.csv")

fig, ax = plt.subplots(figsize=(18, 5))
ax.bar(df["position"], df["mutation_frequency"], width=1.0, color="crimson", alpha=0.7)
ax.set_xlabel("Position (bp)", fontsize=14)
ax.set_ylabel("Mutation Frequency", fontsize=14)
ax.set_title("Mutation Frequency Across blaNDM (K. pneumoniae)", fontsize=16)
ax.set_ylim(0, 1.05)
ax.axhline(y=0.1, color="grey", linestyle="--", alpha=0.5, label="10% threshold")
ax.legend()

known_mutations = {
    "V88L": 262, "M154L": 460, "D130N": 388,
    "A233V": 698, "E152K": 454, "G200D": 599
}
for label, pos in known_mutations.items():
    if pos <= len(df):
        freq = df.loc[df["position"] == pos, "mutation_frequency"].values[0]
        if freq > 0.02:
            ax.annotate(label, (pos, freq), textcoords="offset points",
                       xytext=(0, 10), ha="center", fontsize=8, rotation=45)

plt.tight_layout()
plt.savefig(DATA / "mutation_frequency_barplot.png", dpi=300)
plt.show()
print("Saved: data/mutation_frequency_barplot.png")
