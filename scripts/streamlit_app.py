import streamlit as st
import pandas as pd
import os
import tempfile
import toytree
import toyplot
import toyplot.html
from pathlib import Path

st.set_page_config(layout="wide", page_title="blaNDM Phylogeny Explorer")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

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
    "NDM-97": "#a93226", "NDM-unknown": "#cccccc",
}

def get_variant(tip_label):
    return tip_label.split("_")[0]

def tree_html(tree):
    tmp = tempfile.NamedTemporaryFile(suffix=".html", delete=False)
    tmp.close()
    try:
        toyplot.html.render(tree, tmp.name)
        with open(tmp.name, "r", encoding="utf-8") as f:
            return f.read()
    finally:
        os.unlink(tmp.name)

@st.cache_data
def load_data():
    meta = pd.read_csv(DATA / "sequence_metadata.csv")
    vf = pd.read_csv(RESULTS / "variant_frequencies.csv")
    mf = pd.read_csv(RESULTS / "mutation_frequencies_per_position.csv")
    try:
        am = pd.read_csv(RESULTS / "all_mutations.csv")
    except FileNotFoundError:
        am = pd.DataFrame()
    return meta, vf, mf, am

@st.cache_resource
def load_tree():
    return toytree.tree(DATA / "blaNDM_Kpneumoniae_phylogeny.fasta.contree")

@st.cache_resource
def render_tree_html(_tree):
    labels = _tree.get_tip_labels()
    variants = [get_variant(l) for l in labels]
    colors = [VARIANT_COLORS.get(v, "#cccccc") for v in variants]

    canvas, axes, mark = _tree.draw(
        width=1000,
        height=1800,
        tip_labels_style={"font-size": "10px"},
        tip_labels_colors=colors,
        node_labels=False,
    )
    return tree_html(canvas)

meta, vf, mf, am = load_data()

st.sidebar.title("blaNDM Explorer")
st.sidebar.markdown("Interactive analysis of 250 blaNDM sequences from *K. pneumoniae*")
page = st.sidebar.radio("Page", ["Overview", "Interactive Tree", "Mutation Analysis", "Variant Distribution", "Combined Figure", "Conclusions"])

if page == "Overview":
    st.title("blaNDM Phylogeny & Mutation Analysis")
    st.markdown("Comprehensive analysis of 250 blaNDM sequences from *Klebsiella pneumoniae*")

    total_seqs = 250
    unique_seqs = 90
    alignment_len = 2271
    parsimony_info = 1343
    n_variants = 63
    high_freq_pos = int((mf["mutation_frequency"] >= 0.05).sum())
    top_pos = mf.loc[mf["mutation_frequency"].idxmax()]

    k1, k2, k3, k4 = st.columns(4)
    k1.metric("Total Sequences", total_seqs)
    k2.metric("Unique (Tree)", unique_seqs)
    k3.metric("Alignment Length", f"{alignment_len} bp")
    k4.metric("Parsimony-Informative", parsimony_info)

    k5, k6, k7, k8 = st.columns(4)
    k5.metric("Variant Patterns", n_variants)
    k6.metric("High-Freq Positions (≥5%)", high_freq_pos)
    k7.metric("Most Mutated", f"Pos {int(top_pos['position'])}")
    k8.metric("Top Freq", f"{top_pos['mutation_frequency']:.1%}")

    st.subheader("IQ-TREE Summary")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"""
        - **Model**: HKY+G4 (selected by BIC)
        - **Log-likelihood**: -15,943.42
        - **Bootstrap correlation**: 0.990
        - **Total tree length**: 5.31 substitutions/site
        """)
    with col2:
        st.markdown(f"""
        - **Constant sites**: 529 (23.3%)
        - **Singleton sites**: 399
        - **Near-zero branches (warning)**: 37
        - **Sequences >50% gaps**: 243
        """)

    st.subheader("Key Findings")
    st.markdown(f"""
    1. **Insertion hotspot at positions 789-795**: Present in ~90% of sequences — a key distinguishing feature of most NDM variants vs NDM-1
    2. **Position 1294 (A→C)**: Second most variable at 41.8%, strongly associated with Clade B
    3. **Two major clades**: Clade A (~35.5%, 789-795 insertion) and Clade B (~26.0%, insertion + G941T + A1294C)
    4. **NDM-63 as outgroup**: ~0.68 substitutions/site from the rest — the most divergent sequence in the dataset
    5. **161 identical sequences**: Only 90 unique sequence types in the tree
    6. **Quality concern**: 243 sequences have >50% gaps — likely alignment artifacts from divergent sequences
    """)

    st.subheader("Visualizations")
    c1, c2 = st.columns(2)
    for path, title in [(RESULTS / "mutation_frequency_barplot.png", "Mutation Frequency Bar Plot"),
                         (RESULTS / "mutation_heatmap.png", "Mutation Heatmap")]:
        if path.exists():
            (c1 if title.startswith("Mutation F") else c2).image(str(path), caption=title)

elif page == "Interactive Tree":
    st.title("Interactive Phylogenetic Tree")
    st.markdown("**Consensus tree** from IQ-TREE (90 unique sequence types). Tip colors = NDM variant. Scroll to zoom, drag to pan, hover for labels.")

    with st.spinner("Rendering tree..."):
        tree = load_tree()
        html = render_tree_html(tree)

    st.components.v1.html(html, height=900, scrolling=True)

    with st.expander("Tip Color Legend"):
        cols = st.columns(5)
        for i, (variant, color) in enumerate(sorted(VARIANT_COLORS.items())):
            cols[i % 5].markdown(f'<span style="display:inline-block;width:12px;height:12px;background:{color};border-radius:2px;margin-right:4px"></span> {variant}', unsafe_allow_html=True)

elif page == "Mutation Analysis":
    st.title("Mutation Analysis")

    tab1, tab2, tab3 = st.tabs(["Per-Position Frequency", "Interactive Explorer", "Saved Plots"])

    with tab1:
        st.subheader("Per-Position Mutation Frequency")
        st.markdown(f"**{len(mf)}** positions analyzed. **{(mf['mutation_frequency'] >= 0.05).sum()}** positions with ≥5% mutation frequency.")

        threshold = st.slider("Frequency threshold", 0.0, 1.0, 0.05, 0.01)
        filtered = mf[mf["mutation_frequency"] >= threshold]
        st.write(f"Showing {len(filtered)} positions above {threshold:.0%}")

        chart_data = filtered.rename(columns={"position": "Position", "mutation_frequency": "Frequency"})
        st.area_chart(chart_data.set_index("Position"), height=400)

        st.dataframe(filtered.sort_values("mutation_frequency", ascending=False).head(20),
                     column_config={"position": "Pos", "mutation_frequency": st.column_config.ProgressColumn("Freq", format="%.2f", min_value=0, max_value=1)},
                     use_container_width=True, hide_index=True)

    with tab2:
        st.subheader("Explore by Position Range")
        pos_range = st.slider("Position range", 1, len(mf), (1, 100))
        subset = mf[(mf["position"] >= pos_range[0]) & (mf["position"] <= pos_range[1])]
        st.line_chart(subset.set_index("position"), height=300)

        if not am.empty:
            pos = st.number_input("Show mutations at position", 1, len(mf), 795)
            per_pos = am[am["position"] == pos].head(20)
            if len(per_pos) > 0:
                st.write(f"**Position {pos}**: {per_pos['ref_base'].iloc[0]} → {', '.join(per_pos['query_base'].unique())}")
                st.dataframe(per_pos[["sequence", "query_base"]].head(10), hide_index=True, use_container_width=True)

    with tab3:
        c1, c2 = st.columns(2)
        for path, title in [(RESULTS / "mutation_frequency_barplot.png", "Bar Plot"),
                             (RESULTS / "mutation_heatmap.png", "Heatmap")]:
            if path.exists():
                (c1 if "Bar" in title else c2).image(str(path), caption=title)

elif page == "Variant Distribution":
    st.title("Variant Distribution")

    c1, c2 = st.columns([3, 2])
    with c1:
        top_n = st.slider("Show top N variants", 5, 64, 15)
        chart_vf = vf.sort_values("count", ascending=True).tail(top_n)
        st.bar_chart(chart_vf.set_index("variant")["count"])

    with c2:
        st.metric("Total Variant Patterns", len(vf))
        st.metric("Most Common", vf.iloc[0]["variant"][:30] + "...")
        st.metric("Most Common Count", vf.iloc[0]["count"])
        st.metric("Most Common Frequency", f"{vf.iloc[0]['frequency']:.1%}")

    st.subheader("Variant Table")
    display = vf[["variant", "count", "frequency"]].copy()
    display["variant_short"] = display["variant"].apply(lambda x: x[:80] + "..." if len(x) > 80 else x)
    display["frequency"] = display["frequency"].apply(lambda x: f"{x:.1%}")
    st.dataframe(display[["variant_short", "count", "frequency"]].rename(
        columns={"variant_short": "Variant Mutations", "count": "Count", "frequency": "Freq"}),
        use_container_width=True, hide_index=True)

elif page == "Combined Figure":
    st.title("Combined Figure")
    path = FIGURES / "combined_figure.png"
    if path.exists():
        st.image(str(path), caption="Phylogenetic Tree + Mutation Heatmap", use_container_width=True)
    else:
        st.warning("combined_figure.png not found — run scripts/10_combine_figure.py first")

elif page == "Conclusions":
    st.title("Key Conclusions")

    st.markdown("---")
    c1, c2, c3 = st.columns(3)
    c1.metric("Best-fit Model", "HKY+G4", "by BIC")
    c2.metric("Bootstrap Correlation", "0.990", "1000 replicates")
    c3.metric("Parsimony-Informative", "1,343", "59.1% of sites")

    st.markdown("---")

    with st.container(border=True):
        st.subheader("1. The 789-795 Insertion is Near-Universal")
        st.markdown("""
        **99.2%** of sequences carry a 4-nucleotide insertion (G,G,T,T) at positions 789-795 relative to the NDM-1 reference.
        This insertion occurred early in the divergence of *K. pneumoniae* blaNDM and is now a defining feature of the circulating population —
        only 2 out of 250 sequences match the ancestral NDM-1 pattern.
        """)

    with st.container(border=True):
        st.subheader("2. Two Major Clades Dominate the Population")
        col_a, col_b = st.columns(2)
        with col_a:
            st.markdown("**Clade A** — Insertion Only")
            st.markdown("89 sequences (35.5%)")
            st.markdown("789-795 insertion without additional defining mutations. Includes NDM-5, NDM-7, NDM-50, NDM-29, NDM-9.")
        with col_b:
            st.markdown("**Clade B** — Insertion + G941T + A1294C")
            st.markdown("65 sequences (26.0%)")
            st.markdown("Derived from Clade A with two additional substitutions. Predominantly NDM-1 and NDM-5 sub-variants. Forms a well-supported subcluster.")

    with st.container(border=True):
        st.subheader("3. High Sequence Redundancy")
        st.markdown("""
        **64.4%** of sequences (161/250) are identical after alignment, collapsing to just **90 unique sequence types**.
        57 additional low-frequency variant patterns represent ongoing diversification. This pattern strongly suggests extensive
        clonal expansion of a limited number of highly successful NDM variants rather than independent acquisition events.
        """)

    with st.container(border=True):
        st.subheader("4. NDM-63 is Deeply Divergent")
        st.markdown("""
        The two NDM-63 sequences (OR667647.1, NG_244540.1) form a distinct outgroup with **~0.685 substitutions/site**
        from the main clade — over **100× more divergent** than typical within-clade distances. This suggests either
        an early divergence event predating the main radiation, or a separate acquisition from a distinct bacterial source.
        """)

    with st.container(border=True):
        st.subheader("5. Prototypical NDM-1 is Rare in the Current Population")
        st.markdown("""
        Only **2 out of 250 sequences** (0.8%) match the ancestral NDM-1 pattern (FN396876, first reported in 2009).
        The circulating population has diverged substantially in the intervening years, dominated by the 789-795
        insertion-carrying variants. The reference NDM-1, while useful as a baseline for mutation calling, is not
        representative of contemporary clinical isolates.
        """)

    with st.container(border=True):
        st.subheader("6. Data Quality: Contamination & Limitations")
        col_q1, col_q2 = st.columns(2)
        with col_q1:
            st.markdown("**⚠️ Contamination (~7-10 seqs)**")
            st.markdown("KPC-2, ompK porin genes, and plasmid sequences were misidentified as blaNDM in the NCBI query.")
            st.markdown("**⚠️ >50% gaps: 243 seqs**")
            st.markdown("Expected — indel-rich alignment creates many gap positions relative to NDM-1 reference.")
        with col_q2:
            st.markdown("**⚠️ Composition failure: 19 seqs**")
            st.markdown("Base composition significantly differs from model expectation. May indicate contamination or sequencing bias.")
            st.markdown("**⚠️ No geographic data**")
            st.markdown("Country and year metadata are unavailable. Temporal and phylogeographic analysis is not possible.")

    st.markdown("---")
    st.caption("Report generated from the blaNDM Experiment pipeline | 250 K. pneumoniae genomic DNA sequences | Data analyzed May 2026")
