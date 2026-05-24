import re
import pandas as pd
import toytree
import toyplot
import toyplot.svg
import toyplot.html
import toyplot.pdf

meta = pd.read_csv("data/sequence_metadata.csv", quoting=1)

VARIANT_COLORS = {
    "NDM-1": "#1f77b4", "NDM-4": "#ff7f0e", "NDM-5": "#2ca02c",
    "NDM-7": "#d62728", "NDM-9": "#9467bd", "NDM-10": "#8c564b",
    "NDM-16a": "#e377c2", "NDM-23": "#7f7f7f", "NDM-25": "#bcbd22",
    "NDM-28": "#17becf", "NDM-29": "#aec7e8", "NDM-39": "#ffbb78",
    "NDM-41": "#98df8a", "NDM-43": "#ff9896", "NDM-44": "#c5b0d5",
    "NDM-50": "#c49c94", "NDM-52": "#f7b6d2", "NDM-54": "#c7c7c7",
    "NDM-61": "#dbdb8d", "NDM-63": "#9edae5", "NDM-71": "#f5b041",
    "NDM-72": "#5dade2", "NDM-73": "#e74c3c", "NDM-84": "#27ae60",
    "NDM-85": "#8e44ad", "NDM-94": "#d35400", "NDM-96": "#2980b9",
    "NDM-97": "#a93226",
}

def get_variant(seq_id):
    if "FN396876" in seq_id:
        return "NDM-1"
    accession = seq_id.split("_", 1)[1] if "_" in seq_id else seq_id
    row = meta[meta["accession"] == accession]
    if len(row) == 0:
        return "NDM-unknown"
    original = row.iloc[0]["original_id"]
    m = re.search(r"NDM-?\d+", original)
    return m.group(0) if m else "NDM-unknown"

tree = toytree.tree("data/blaNDM_Kpneumoniae_phylogeny.fasta.contree")

tip_colors = []
for tip in tree.get_tip_labels():
    variant = get_variant(tip)
    tip_colors.append(VARIANT_COLORS.get(variant, "#cccccc"))

canvas, axes, mark = tree.draw(
    width=1200,
    height=1800,
    tip_labels=True,
    tip_labels_colors=tip_colors,
    node_labels="support",
    node_labels_style={"font-size": "8px"},
    node_sizes=8,
)

toyplot.svg.render(canvas, "figures/phylogeny_tree.svg")
toyplot.html.render(canvas, "figures/phylogeny_tree.html")
toyplot.pdf.render(canvas, "figures/phylogeny_tree.pdf")
print("Saved figures/phylogeny_tree.{svg,html,pdf}")
