from Bio import SearchIO
import os
def parse_query_results(path):
    """
    This fucntion will p[arse the query results
    :param path: The path of the file that contains the query results
    :return: An array that contains the results described in the file
    """
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

def get_pairs(query_results):
    """
    This function will build and return the paisrs of (miRNA,host) for every full match
    :param query_results:The given array of results
    :return: An array of pairs
    """
    pairs = []
    for res in query_results:
        pident = float(res["pident"])
        if pident == 100:
            miRNA = res['qseqid']
            host = res["sseqid"]

            try:
                index_semicolon = host.index(";")
                index = host.index("|")
                pre = host[:index_semicolon]
                post = host[index]
                host = "%s%s" % (pre,post)

            except(ValueError):
                pass
            pairs.append((miRNA, host))
    return pairs

def get_3_tag_UTR(host_header):
    """
    This function will return the sequence that the given header describes (in the mRNA file)
    :param host_header: The given host header
    :return: The sequence that the given header describes (in the mRNA file)
    """
    mRNA_path_file = "%s\%s" %(os.path.dirname(os.path.abspath(__file__)),r"Resources\mRNA.fasta")
    mRNA_file = open(mRNA_path_file,"r")
    string_to_find = ">%s" % host_header
    start_sequence_assembly = False
    stop_search = False
    sequence = ""

    for line in mRNA_file:
        if stop_search:
            return sequence
        if start_sequence_assembly:
                if line[0] == ">":
                    stop_search = True
                else:
                    sequence = "%s%s" % (sequence,line[:-1])
        else:
            if string_to_find in line:
                start_sequence_assembly = True

    return None # Should not get here

def get_host_sequences(pairs):
    """
    This function will find the sequences of all the host mRNAs and will add them to the pairs array
    :param pairs: [(miRNA,mRNA)]
    :return: [(miRNA,(mRNA,sequence))]
    """
    new_pairs = []
    for pair in pairs:
        host_header = pair[1]
        host_sequence = get_3_tag_UTR(host_header)
        new_pairs.append((pair[0],(host_header,host_sequence)))
    return new_pairs

path = "%s\%s" %(os.path.dirname(os.path.abspath(__file__)),r"Resources\res_blastn_compact.txt")
query_results = parse_query_results(path)
pairs = get_pairs(query_results)
print(len(pairs))

n_pairs = get_host_sequences(pairs)
for pair in n_pairs:
    print(pair)
"""
res = get_3_tag_UTR("WBGene00001520|WBGene00001520.1|K09A9.5.1|gas-1")
print(res)
seq = "ATACTCAACTCATTAGGCACGTAGACGGATTCTCTATAGCACATTCTCAACTCACTCTTTATTTATCCCTTGCACCGCGAAATGTGTTCTGTTTATTTTCTTTTCTTTTCAAAGTTCGTTTTTTTCTGAAATTCAAAAATGTAGATTTGTATCTGTTCACTGTTGTTCAGGATGTTCAATAAATAAGTCTGCAACCAAGACGCAAAG"
print(seq == res)
"""


