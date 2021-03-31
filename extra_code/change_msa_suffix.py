from Bio import SeqIO

msa_path = "/Users/noa/Workspace/data/LARGE_FILES/information_500gene_Pro.nex"
destination_path = "/Users/noa/Workspace/data/LARGE_FILES/information_500gene_Pro.fasta"

with open(msa_path) as input_handle, open(
   destination_path, "w"
) as output_handle:
    sequences = SeqIO.parse(input_handle, "nexus")
    count = SeqIO.write(sequences, output_handle, "fasta")

print("Converted %i records" % count)