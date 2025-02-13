# /Users/guillermocomesanacimadevila/raw_genome_test/ncbi_dataset/data
# gene.fna

import os

def read_file(location):
    with open(os.path.expanduser(location), 'r') as file:
        return file.read()

def count_nuc(seq):
    return print(f"A: {seq.count("A")}, T: {seq.count("T")}, C: {seq.count("C")}, G: {seq.count("G")}, U: {seq.count("U")}")

def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)

gene = read_file("/Users/guillermocomesanacimadevila/raw_genome_test/ncbi_dataset/data/gene.fna")
indices = [i for i, x in enumerate(gene) if x == "U"]

print(len(gene))
print(count_nuc(gene))
print(gc_content(gene), f"G & C: {gc_content(gene) / 2}")
print(indices) # 2Us # fine not shitty Us -> part of the analysis