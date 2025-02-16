import os
import matplotlib.pyplot as plt

def read_file(location):
    with open(os.path.expanduser(location), "r") as file:
        sequence = ""
        for line in file:
            line = line.strip()
            if not line.startswith(">"):
                sequence += line
    return sequence

# take each position == 50bp long
# compute gc_content
# accumulate start position

def gc_content(seq, kmer):
    gc = []
    for pos in range(len(seq) - kmer + 1):
        sub_seq = seq[pos:pos + kmer]
        gc_count = sub_seq.count("G") + sub_seq.count("C")
        gc.append((gc_count / kmer) * 100)
    return gc

duox2 = read_file("/Users/guillermocomesanacimadevila/Desktop/gene_testing/ncbi_dataset/data/gene.fna")
gc_cont = {position: gc for position, gc in enumerate(gc_content(duox2, 50))}
print(gc_cont)

# plot
plt.figure(figsize=(10, 6))

plt.plot([key for key, value in gc_cont.items()], [val for key, val in gc_cont.items()], marker="o", linestyle="-", linewidth=2, color="#228B22")
plt.xlabel("Position", fontsize=14, labelpad=10)
plt.ylabel("GC Content (%)", fontsize=14, labelpad=10)
plt.title("DUOX2-wide GC Content", fontsize=16, pad=15)
plt.legend(fontsize=12, loc="upper left", frameon=True, shadow=True)

plt.yticks(fontsize=12)
plt.grid(visible=True, linestyle="--", alpha=0.7)

plt.tight_layout()
plt.show()