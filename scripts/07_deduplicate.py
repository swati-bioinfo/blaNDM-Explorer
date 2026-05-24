from Bio import SeqIO

records = list(SeqIO.parse("blaNDM_Kpneumoniae_aligned_trimmed.fasta", "fasta"))
seen_seqs = {}
unique = []
for rec in records:
    seq_str = str(rec.seq).upper()
    if seq_str not in seen_seqs:
        seen_seqs[seq_str] = rec.id
        unique.append(rec)

SeqIO.write(unique, "blaNDM_Kpneumoniae_aligned_dedup.fasta", "fasta")
print(f"Deduplicated: {len(records)} → {len(unique)} sequences")