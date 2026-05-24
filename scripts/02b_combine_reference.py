from Bio import SeqIO
from pathlib import Path

DATA = Path("data")

aln = list(SeqIO.parse(DATA / "blaNDM_Kpneumoniae_aligned.fasta", "fasta"))
ref = list(SeqIO.parse(DATA / "blaNDM_NDM1_reference.fasta", "fasta"))

ref[0].id = "NDM-1_FN396876"
ref[0].description = ""

aln.append(ref[0])
SeqIO.write(aln, DATA / "blaNDM_Kpneumoniae_wref.fasta", "fasta")
print(f"Combined: {len(aln)} sequences (original {len(aln)-1} + NDM-1 reference) in data/blaNDM_Kpneumoniae_wref.fasta")