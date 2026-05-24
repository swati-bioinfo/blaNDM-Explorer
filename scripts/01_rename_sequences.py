#!/usr/bin/env python3
import csv
from Bio import SeqIO
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"

input_file = DATA / "blaNDM_Kpneumoniae_raw.fasta"
output_file = DATA / "blaNDM_Kpneumoniae_renamed.fasta"
metadata_file = DATA / "sequence_metadata.csv"

records = []
rows = []
for i, record in enumerate(SeqIO.parse(input_file, "fasta")):
    header = record.description
    accession = header.split()[0] if " " in header else header

    parts = header.split("|")
    organism = parts[1] if len(parts) > 1 else "Unknown"
    country = parts[3] if len(parts) > 3 else "Unknown"
    year = parts[4] if len(parts) > 4 else "Unknown"

    short_id = f"KpNDM_{i+1:04d}_{country}_{year}"
    record.id = short_id
    record.description = ""
    records.append(record)
    rows.append([header, short_id, accession, organism, country, year])

with open(metadata_file, "w", newline="") as f:
    writer = csv.writer(f, quoting=csv.QUOTE_ALL)
    writer.writerow(["original_id", "short_id", "accession", "organism", "country", "year"])
    writer.writerows(rows)

SeqIO.write(records, output_file, "fasta")
print(f"Renamed {len(records)} sequences -> {output_file}")