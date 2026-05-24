import re
import pandas as pd
from Bio import SeqIO

meta = pd.read_csv("data/sequence_metadata.csv")
id_map = dict(zip(meta["short_id"], meta["original_id"]))

records = []
for rec in SeqIO.parse("data/blaNDM_Kpneumoniae_aligned_trimmed.fasta", "fasta"):
    if "FN396876" in rec.id:
        # NDM-1 reference — keep as-is (already NDM-1_FN396876)
        new_id = rec.id
    else:
        original_header = id_map.get(rec.id, rec.id)
        accession = original_header.split()[0] if " " in original_header else original_header
        match = re.search(r"NDM-\d+", original_header)
        ndm_variant = match.group(0) if match else "NDM-unknown"
        new_id = f"{ndm_variant}_{accession}"
    rec.id = new_id
    rec.description = ""
    records.append(rec)

SeqIO.write(records, "data/blaNDM_Kpneumoniae_phylogeny.fasta", "fasta")
print(f"Renamed {len(records)} sequences -> data/blaNDM_Kpneumoniae_phylogeny.fasta")