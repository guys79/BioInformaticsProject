from Bio import SearchIO


def parse_query_results(path):
    file = open(path,"r")
    results = []
    for line in file:
        split_by_tab = line.split("\t")
        dictionary = {}
        dictionary["qseqid"] = split_by_tab[0]
        dictionary["sseqid"] = split_by_tab[1]
        dictionary["pident"] = split_by_tab[2]
        dictionary["length"] = split_by_tab[3]
        dictionary["mismatch"] = split_by_tab[4]
        dictionary["gapopen"] = split_by_tab[5]
        dictionary["qstart"] = split_by_tab[6]
        dictionary["qend"] = split_by_tab[7]
        dictionary["sstart"] = split_by_tab[8]
        dictionary["send"] = split_by_tab[9]
        dictionary["evalue"] = split_by_tab[10]
        dictionary["bitscore"] = split_by_tab[11][:-1]

        results.append(dictionary)
    return results

path = r"C:\Users\guys79\Downloads\res_blastn_compact.txt"
query_results = parse_query_results(path)
for res in query_results:
    print(res)
