from Bio import SeqIO

# TO check the path of vector and gene of interest(GOI)
vector_path = "fasta_files/pCAMBIA1301.fasta"
goi_path = "fasta_files/RAS1.fasta"

# Reading FASTA sequence of vector
vector_record = SeqIO.read(vector_path, "fasta")
print(f"Vector ID: {vector_record.id}")
print(f"Vector Length: {len(vector_record.seq)}")

# Read FASTA sequence gene of interest (GOI)
goi_record = SeqIO.read(goi_path, "fasta")
print(f"GOI ID: {goi_record.id}")
print(f"GOI Length: {len(goi_record.seq)}")
# Define an insertion position
insertion_site = 3000
modified_vector = vector_record.seq[:insertion_site] + goi_record.seq + vector_record.seq[insertion_site:]

# Modifying FASTA sequence
from Bio.SeqRecord import SeqRecord
modified_record = SeqRecord(modified_vector, id="Modified_pCAMBIA1301", description="pCAMBIA1301 with GOI Insert")

# Saving the modified sequences
output_path = "fasta_files/Modified_pCAMBIA1301.fasta"
SeqIO.write(modified_record, output_path, "fasta")
print("Modified plasmid saved successfully.")
from Bio.Restriction import *
from Bio.Seq import Seq

# Convert sequences to Biopython Seq objects
vector_seq = vector_record.seq
goi_seq = goi_record.seq

# Restriction enzymes
enzymes = [enzyme for enzyme in AllEnzymes]
# Check cut sites in vector
print("Restriction Sites in Vector:")
for enzyme in enzymes:
    cut_sites = enzyme.search(vector_seq)
    if cut_sites:
        print(f"{enzyme.__name__} cuts at positions: {cut_sites}")

# Checking the for sites in GOI
print("\nRestriction Sites in GOI:")
for enzyme in enzymes:
    cut_sites = enzyme.search(goi_seq)
    if cut_sites:
        print(f"{enzyme.__name__} cuts at positions: {cut_sites}")
# EcoRI cut site in vector (assume at position 1000)
vector_cut_pos = 1000
vector_fragments = [vector_seq[:vector_cut_pos], vector_seq[vector_cut_pos:]]

# Simulate EcoRI digestion of GOI
goi_cut_pos = 1000
goi_fragments = [goi_seq[:goi_cut_pos], goi_seq[goi_cut_pos:]]

print("Vector and GOI digested successfully.")
# Simulate ligation by inserting GOI into vector at cut position
recombinant_plasmid = vector_fragments[0] + goi_fragments[1] + vector_fragments[1]

print(f"Recombinant Plasmid Length: {len(recombinant_plasmid)} bp")
from Bio.SeqRecord import SeqRecord

# Create a new sequence record
recombinant_record = SeqRecord(recombinant_plasmid, id="Recombinant_pCAMBIA1301", description="pCAMBIA1301 with GOI")

# Save as FASTA
output_path = "fasta_files/Recombinant_pCAMBIA1301.fasta"
SeqIO.write(recombinant_record, output_path, "fasta")

print("Recombinant plasmid saved successfully.")

import dna_features_viewer as dfv
from Bio import SeqIO
import matplotlib.pyplot as plt

# Load recombinant plasmid sequence
recombinant_path = "fasta_files/Recombinant_pCAMBIA1301.fasta"
recombinant_record = SeqIO.read(recombinant_path, "fasta")

# Define plasmid length
plasmid_length = len(recombinant_record.seq)

# GOI insertion site
goi_start = 6447  # Restriction enzyme cut site
goi_end = goi_start + 2904  # GOI length

# Define features (plasmid backbone, GOI, restriction site)
plasmid_backbone = dfv.GraphicFeature(start=0, end=plasmid_length, strand=0, color="blue", label="Plasmid Backbone")
goi_feature = dfv.GraphicFeature(start=goi_start, end=goi_end, strand=+1, color="red", label="GOI Insert")
restriction_site = dfv.GraphicFeature(start=6447, end=6457, strand=0, color="green", label="BstZ17I Site")

# Create plasmid map
diagram = dfv.CircularGraphicRecord(sequence_length=plasmid_length, features=[plasmid_backbone, goi_feature, restriction_site])

# Plot circular plasmid
fig, ax = plt.subplots(figsize=(6, 6))
diagram.plot(ax=ax, figure_width=10)
plt.show()
