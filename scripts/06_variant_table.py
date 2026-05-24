import pandas as pd
from Bio import SeqIO
from pathlib import Path

DATA = Path("data")

aln = list(SeqIO.parse(DATA / "blaNDM_Kpneumoniae_aligned_trimmed.fasta", "fasta"))

ref_rec = None
for rec in aln:
    if "FN396876" in rec.id or "NDM_0001" in rec.id or "NDM-1" in rec.id:
        ref_rec = rec
        break

if ref_rec is None:
    ref_rec = aln[0]
    print(f"Using {ref_rec.id} as reference (NDM-1 not found; using first sequence)")

ref_seq = str(ref_rec.seq).upper()

variants = {}
for rec in aln:
    seq = str(rec.seq).upper()
    changes = []
    for i, (r, q) in enumerate(zip(ref_seq, seq)):
        if r != q and q not in ["-", "N"]:
            changes.append(f"{r}{i+1}{q}")
    if not changes:
        variant_key = "Reference (NDM-1 / most common)"
    else:
        variant_key = "+".join(changes)

    if variant_key not in variants:
        variants[variant_key] = {"count": 0, "sequences": []}
    variants[variant_key]["count"] += 1
    variants[variant_key]["sequences"].append(rec.id)

df_variants = pd.DataFrame([
    {"variant": k, "count": v["count"], "frequency": v["count"]/len(aln),
     "example_sequences": "; ".join(v["sequences"][:3])}
    for k, v in sorted(variants.items(), key=lambda x: -x[1]["count"])
])

df_variants.to_csv(DATA / "variant_frequencies.csv", index=False)
print(df_variants.head(15).to_string())
print(f"\nTotal unique variant patterns: {len(df_variants)}")
print(f"Saved: {DATA / 'variant_frequencies.csv'}")
