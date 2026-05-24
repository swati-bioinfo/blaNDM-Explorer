from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

aln = list(SeqIO.parse("data/blaNDM_Kpneumoniae_aligned_trimmed.fasta", "fasta"))

ref_rec = None
for rec in aln:
    if "FN396876" in rec.id or "NDM_0001" in rec.id or "NDM-1" in rec.id:
        ref_rec = rec
        break

if ref_rec is None:
    ref_rec = aln[0]
    print(f"Using {ref_rec.id} as reference (NDM-1 not found; using first sequence)")

ref_seq = str(ref_rec.seq).upper()
aln_length = len(ref_seq)

mutations = []
mutation_matrix = []

for rec in aln:
    seq = str(rec.seq).upper()
    row = []
    seq_muts = []
    for pos in range(aln_length):
        ref_base = ref_seq[pos]
        query_base = seq[pos]
        if ref_base != query_base and query_base not in ["-", "N"]:
            row.append(1)
            seq_muts.append({
                "sequence": rec.id,
                "position": pos + 1,
                "ref_base": ref_base,
                "query_base": query_base,
                "is_mutation": True
            })
        else:
            row.append(0)
    mutation_matrix.append(row)
    mutations.extend(seq_muts)

df_mutations = pd.DataFrame(mutations)
df_mutation_matrix = pd.DataFrame(mutation_matrix)

pos_frequencies = df_mutation_matrix.mean(axis=0)
df_pos_freq = pd.DataFrame({
    "position": range(1, aln_length + 1),
    "mutation_frequency": pos_frequencies.values
})

df_pos_freq.to_csv("data/mutation_frequencies_per_position.csv", index=False)
df_mutations.to_csv("data/all_mutations.csv", index=False)

print(f"Total mutations detected: {len(df_mutations)}")
print(f"Variable positions: {(pos_frequencies > 0).sum()}/{aln_length}")
print(f"Most variable position: Position {df_pos_freq.loc[df_pos_freq['mutation_frequency'].idxmax(), 'position']} "
      f"(freq = {df_pos_freq['mutation_frequency'].max():.3f})")