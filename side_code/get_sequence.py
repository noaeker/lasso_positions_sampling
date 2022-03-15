from Bio import SeqIO

original_alignment_path = "/Users/noa/Workspace/data/ABC_DR/Selectome/Euteleostomi/ENSGT00680000099951/example_msa.phy"
with open(original_alignment_path) as original:
    original_alignment_data = list(SeqIO.parse(original, 'phylip-relaxed'))

a = str(original_alignment_data[5].seq)
print(a.replace('-',''))

print(type(original_alignment_data[0].seq))
