from Bio import SeqIO

msa_path = "/Users/noa/Workspace/data/LARGE_FILES/Shen_et_al_2016/308AA_combined.fasta"
destination_path = "/Users/noa/Workspace/data/LARGE_FILES/information_500gene_Pro.fasta"


records = list(SeqIO.parse(msa_path, "fasta"))
print(len(records))
print(len(records[0].seq))
# with open(msa_path) as input_handle, open(
#    destination_path, "w"
# ) as output_handle:
#     sequences = SeqIO.parse(input_handle, "nexus")
#     count = SeqIO.write(sequences, output_handle, "fasta")
#
# print("Converted %i records" % count)