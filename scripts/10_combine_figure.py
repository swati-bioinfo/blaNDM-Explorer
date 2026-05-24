import os
import fitz
from PIL import Image

if not os.path.exists("figures/phylogeny_tree.pdf"):
    print("No phylogeny_tree.pdf found — run 09_tree_visualize.py first")
    exit(1)

tree_doc = fitz.open("figures/phylogeny_tree.pdf")
tree_pix = tree_doc[0].get_pixmap(dpi=200)
tree_pix.save("figures/phylogeny_tree.png")
tree_doc.close()

heatmap_path = "figures/mutation_heatmap.png"
if not os.path.exists(heatmap_path):
    heatmap_path = "results/mutation_heatmap.png"

if not os.path.exists(heatmap_path):
    print("No mutation_heatmap.png found (checked figures/ and results/)")
    print("Run 05_plot_mutation_heatmap.py to generate it")
    exit(0)

tree_img = Image.open("figures/phylogeny_tree.png")
heatmap_img = Image.open(heatmap_path)

tree_w, tree_h = tree_img.size
heatmap_w, heatmap_h = heatmap_img.size
new_heatmap_h = tree_h
new_heatmap_w = int(heatmap_w * new_heatmap_h / heatmap_h)
heatmap_img = heatmap_img.resize((new_heatmap_w, new_heatmap_h), Image.LANCZOS)

combined = Image.new("RGB", (tree_w + new_heatmap_w, tree_h))
combined.paste(tree_img, (0, 0))
combined.paste(heatmap_img, (tree_w, 0))
combined.save("figures/combined_figure.png", dpi=(300, 300))
print("Saved figures/combined_figure.png")
