from Bio.Seq import Seq
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
    This function will build and return the pairs of (miRNA,host) for every full match
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

def get_3_tag_UTR(host_header,file_name):
    """
    This function will return the sequence that the given header describes (in the mRNA file)
    :param host_header: The given host header
    :return: The sequence that the given header describes (in the mRNA file)
    """
    mRNA_path_file = r"%s\Resources\%s.fasta" %(os.path.dirname(os.path.abspath(__file__)),file_name)
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
    if sequence!= "":
        return sequence

    return None # Should not get here

def get_cbr_miRNA_That_Met_threshold(query_results_cbr):
    """
    This function will add the 'is_conserved' key to the dictionaries
    :param query_results_cbr: the list of dictionaries
    :return: None
    """
    for query_res in query_results_cbr:
        pident = float(query_res['pident'])
        if  pident == 100:
            query_res['is_conserved'] = True
        else:
            query_res['is_conserved'] = False

def get_miRNA_and_host_sequences(pairs,miRNA_flie_name,mRNA_file_name):
    """
    This function will find the sequences of all the host mRNAs and will add them to the pairs array
    :param pairs: [(miRNA,mRNA)]
    :return: [(miRNA,(mRNA,sequence))]
    """
    new_pairs = []
    leng = len(pairs)
    count = 1
    for pair in pairs:
        print("%d/%d" % (leng,count))
        count+=1
        host_header = pair[1]
        host_sequence = get_3_tag_UTR(host_header,mRNA_file_name)
        miRNA_header = pair[0]
        miRNA_sequence = get_3_tag_UTR(miRNA_header, miRNA_flie_name)
        new_pairs.append(((miRNA_header,miRNA_sequence),(host_header,host_sequence)))
    return new_pairs


def get_seed(sequence):
    """
    This function will return the seed of the sequence
    :param sequence: The sequence
    :return: The seed of the sequence
    """
    return sequence[1:8]


def is_seed_criteria_met(miRNA,mRNA):
    """
    This function will check if the seed criteria is met
    :param miRNA: The miRNA sequence
    :param mRNA: The mRNA sequence
    :return: True IFF the seed criteria is met
    """
    if mRNA == "Sequence unavailable":
        return "unknown"
    miRNA_seed = get_seed(miRNA)
    miRNA_seed = miRNA_seed.replace("U","T")
    miRNA_seq = Seq(miRNA_seed)
    miRNA_seed = miRNA_seq.reverse_complement()
    return miRNA_seed._data in mRNA

def is_seed_criteria_met_in_pairs(pairs):
    """
    This function will return a dictionary that contains weather the seed criteria is met.
    This process occurs for every pair in the pairs' list
    :param pairs: The given pairs list
    :return: A dictionary that contains weather the seed criteria is met.
    """
    res = []
    for pair in pairs:
        dictionary = {}
        miRNA_seq = pair[0][1]
        mRNA_seq = pair[1][1]
        is_target = is_seed_criteria_met(miRNA_seq,mRNA_seq)
        dictionary["miRNA_name"] = pair[0][0]
        dictionary["miRNA_seq"] = pair[0][1]
        dictionary["Host_name"] = pair[1][0]
        dictionary["Host_seq"] = pair[1][1]
        dictionary["is_target"] = is_target
        res.append(dictionary)
    return res


def parse_results_to_csv(query_results):
    """
    This function will parse the results and save them in a csv file name 'results.csv'
    :param query_results: The given data to save
    :return: None
    """
    with open('results.csv', 'w') as file:
        file.write("")
    with open('results.csv','a') as file:
        header = ""
        for key in query_results[0].keys():
            header = "%s%s," % (header,key)
        header = "%s\n" % header
        file.write(header)

    for parse_dictionary in query_results:
        with open('results.csv','a') as f:
            line = ""
            for key in parse_dictionary.keys():
                line = "%s%s," % (line, parse_dictionary[key])
            line = "%s\n" % line
            f.write(line)


def get_all_cell_dictionaries(path_cel,path_cel_pre):
    """
    This function will return all the C.elegans miRNA's
    :param path_cel: The path to the C.elegans miRNA's files
    :param path_cel_pre: The path to the pre-mature C.elegans miRNA's files
    :return: A dictionary containing all the C.elegans miRNA's
    """
    cel_file = open(path_cel,'r')

    all_cell_dictionary = []
    odd = True
    name = None
    seq = None
    for line in cel_file:
        if odd:
            name = line[1:-1]
            odd = False
        else:
            cell_dictionary = {}
            seq = line[:-1]
            odd = True
            cell_dictionary['C.elegans mature name'] = name
            cell_dictionary['C.elegans mature sequence'] = seq
            all_cell_dictionary.append(cell_dictionary)

    odd = True
    name = None
    seq = None
    cel_file_pre = open(path_cel_pre, 'r')
    all_cell_pre_dictionary = []
    for line in cel_file_pre:
        if odd:
            name = line[1:-1]
            odd = False
        else:
            cell_dictionary = {}
            seq = line[:-1]
            odd = True
            cell_dictionary['C.elegans pre-miRNA name'] = name
            cell_dictionary['C.elegans pre-mRNA sequence'] = seq
            all_cell_pre_dictionary.append(cell_dictionary)

    combined = []

    for i in range(len(all_cell_dictionary)):
        for j in range(len(all_cell_pre_dictionary)):
            pre_name = all_cell_pre_dictionary[j]['C.elegans pre-miRNA name']
            pre_seq = all_cell_pre_dictionary[j]['C.elegans pre-mRNA sequence']
            cel_name = all_cell_dictionary[i]['C.elegans mature name']
            cel_seq = all_cell_dictionary[i]['C.elegans mature sequence']
            if pre_name[:pre_name.rindex('_')] == cel_name[:cel_name.rindex('_')]:
                dict = {}
                dict['C.elegans pre-miRNA name'] = pre_name
                dict['C.elegans pre-miRNA sequence'] = pre_seq
                dict['C.elegans mature name'] = cel_name
                dict['C.elegans mature sequence'] = cel_seq
                combined.append(dict)

    return combined

def add_host_data(all_cell_dictionary, final_pairs_mRNA):
    """
    This function will add the host data to the given cel dictionary
    :param all_cell_dictionary: The cell dictionary
    :param final_pairs_mRNA: The host data
    :return: None
    """
    for i in range(len(all_cell_dictionary)):
        cell_dictionary = all_cell_dictionary[i]
        name = cell_dictionary['C.elegans mature name']
        host_exists = False
        for j in range(len(final_pairs_mRNA)):
            if final_pairs_mRNA[j]["miRNA_name"] == name:
                host_exists = True
                host_name = final_pairs_mRNA[j]["Host_name"]
                is_target = final_pairs_mRNA[j]["is_target"]
                cell_dictionary["Host gene name"] = host_name
                if is_target == True:
                    cell_dictionary["Targets the host gene"] = 'yes'
                elif is_target == False:
                    cell_dictionary["Targets the host gene"] = 'no'
                else:
                    cell_dictionary["Targets the host gene"] = 'unknown'
        if not host_exists:
            cell_dictionary["Host gene name"] = "-"
            cell_dictionary["Targets the host gene"] = '-'


def add_cbr_data(all_cell_dictionary, query_results_cbr):
    """
    This function will add the C.briggsae data to the given cel dictionary
    :param all_cell_dictionary: The cell dictionary
    :param query_results_cbr: The C.briggsae data
    :return: None
    """
    for i in range(len(all_cell_dictionary)):
        cell_dictionary = all_cell_dictionary[i]
        name = cell_dictionary['C.elegans mature name']
        cbr_exists = False
        for j in range(len(query_results_cbr)):
            if query_results_cbr[j]["qseqid"] == name:
                cbr_exists = True
                cbr_name = query_results_cbr[j]["sseqid"]
                if query_results_cbr[j]["is_conserved"]:
                    cell_dictionary["Conserved in C.briggsae"] = cbr_name
                else:
                    cell_dictionary["Conserved in C.briggsae"] = False
        if not cbr_exists:
            cell_dictionary["Conserved in C.briggsae"] = "-"



if __name__ == "__main__":

    path_mRNA = "%s\%s" %(os.path.dirname(os.path.abspath(__file__)),r"Resources\res_blastn_compact_mRNA.fasta")
    path_cbr = "%s\%s" %(os.path.dirname(os.path.abspath(__file__)),r"Resources\res_blastn_compact_cbr.fasta")
    path_cel = "%s\%s" %(os.path.dirname(os.path.abspath(__file__)),r"Resources\cel.fasta")
    path_cel_pre = "%s\%s" %(os.path.dirname(os.path.abspath(__file__)),r"Resources\cel-pre.fasta")

    print("parsing query result")
    query_results_mRNA = parse_query_results(path_mRNA)
    query_results_cbr = parse_query_results(path_cbr)

    print("get pairs")
    pairs_mRNA = get_pairs(query_results_mRNA)
    get_cbr_miRNA_That_Met_threshold(query_results_cbr)

    print("get host")
    new_pairs_mRNA = get_miRNA_and_host_sequences(pairs_mRNA,"cel","mRNA")

    print("Updating seed criteria")
    final_pairs_mRNA = is_seed_criteria_met_in_pairs(new_pairs_mRNA)

    print("Gathering data")
    all_cell_dictionary = get_all_cell_dictionaries(path_cel,path_cel_pre)
    add_host_data(all_cell_dictionary,final_pairs_mRNA)
    add_cbr_data(all_cell_dictionary,query_results_cbr)
    parse_results_to_csv(all_cell_dictionary)