from Bio import SeqIO
import numpy as np
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
aln = list(SeqIO.parse(DATA / "blaNDM_Kpneumoniae_aligned_trimmed.fasta", "fasta"))
print(f"Sequences: {len(aln)}")
print(f"Alignment length: {len(aln[0].seq)}")

gap_counts = [str(rec.seq).count("-") for rec in aln]
print(f"Gaps per sequence — Mean: {np.mean(gap_counts):.1f}, Max: {max(gap_counts)}")
print(f"Sequences with >50% gaps: {sum(1 for g in gap_counts if g/len(aln[0].seq) > 0.5)}")

from collections import Counter
conserved = 0
for i in range(len(aln[0].seq)):
    col = [str(rec.seq)[i] for rec in aln]
    if col.count(col[0]) == len(col):
        conserved += 1
print(f"Fully conserved positions: {conserved}/{len(aln[0].seq)} ({100*conserved/len(aln[0].seq):.1f}%)")