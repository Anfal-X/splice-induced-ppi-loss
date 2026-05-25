import pandas as pd
from Bio import SeqIO
Final_table=pd.read_csv('combined_ensemble_input_H.csv')
print("Columns in your table:", Final_table.columns.tolist())
print("Total rows:", len(Final_table))
unique_genes = Final_table["Gene_Symbol"].unique().tolist()
print(f"Number of unique genes: {len(unique_genes)}")
print("First 10 genes:", unique_genes[:10])
P_S='protein_sequences.fasta'
all_sequences={} # holds all sequences
for record in SeqIO.parse(P_S,'fasta'):
    if record.id.endswith("_1"):           # only keep canonical isoform
        gene_name = record.id.rsplit("_", 1)[0]   # strips "_1" → gives "ABCB7"
    sequence=str(record.seq) # record.seq is a special BioPython object representing the sequence.
    all_sequences[gene_name]=sequence
print(f"Sequences loaded from FASTA: {len(all_sequences)}")
filtered_sequences = {} # holds sequences of the genes that are present in final table
missing_genes = []
for gene in unique_genes:
    if gene in all_sequences:
        filtered_sequences[gene] = all_sequences[gene]
    else:
        missing_genes.append(gene)
print(f"Sequences matched: {len(filtered_sequences)}")
print(f"Genes with no sequence in FASTA: {len(missing_genes)}")
if missing_genes:
    print("Missing genes:", missing_genes)
output_fasta = "canonical_sequences_filtered.fasta"
with open(output_fasta, "w") as f:
    for gene, sequence in filtered_sequences.items():
        f.write(f">{gene}\n")
        f.write(f"{sequence}\n")
print(f"Done! Output written to: {output_fasta}")


