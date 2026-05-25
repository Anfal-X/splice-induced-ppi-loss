from Bio import SeqIO

input_file = "ORFs_sequences.fasta"
output_file = "protein_sequences.fasta"

with open(output_file, "w") as out_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        protein = record.seq.translate(to_stop=True)
        out_handle.write(f">{record.id}\n{protein}\n")
