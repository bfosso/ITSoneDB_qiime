import os
from string import strip
import argparse
# import argcomplete
import subprocess
import shlex


def options():
    parser = argparse.ArgumentParser(description="Format ITSoneDB fasta sequence and taxonomy for QIIME",
                                     prefix_chars="-")
    parser.add_argument("-f", "--itsonedb_total_fasta", type=str,
                        help="file containing all the ITSoneDB fasta sequences",
                        action="store", required=True)
    parser.add_argument("-t", "--taxa_folder", type=str,
                        help="folder containing the NCBI taxonomy data",
                        action="store", required=True)
    parser.add_argument("-p", "--identity_clustering", type=float,
                        help="Identity percent for clustering [default = 0.99]",
                        action="store", required=False, default=0.99)
    # argcomplete.autocomplete(parser)
    return parser.parse_args()



def clustering_execution(fasta, id_perc):
    cluster = "vsearch_cluster_%i" % int(id_perc * 100)
    centroids = "%i_ref_ITSoneDB.fa" % int(id_perc * 100)
    cmd = shlex.split("vsearch --cluster_fast %s -uc %s --centroids %s --id %f" % (fasta, cluster, centroids, id_perc))
    p = subprocess.Popen(cmd)
    p.wait()
    cluster_dict = {}
    with open(cluster) as a:
        for line in a:
            s = map(strip, line.split("\t"))
            if s[0] == "S":
                cluster_dict.setdefault(acc_name(s[-2]), [])
                cluster_dict[acc_name(s[-2])].append(acc_name(s[-2]))
            elif s[0] == "H":
                cluster_dict[acc_name(s[-1])].append(acc_name(s[-2]))
    print len(cluster_dict)
    return cluster_dict, centroids


def acc_name(l):
    return l.split("|")[0]


def generate_taxonomic_info(taxa_folder, itsonedb_fasta):
    acc2node = {}
    with open(itsonedb_fasta) as fasta:
        for line in fasta:
            if line.startswith(">"):
                s = map(strip, line.split("|"))
                acc2node[s[0].strip(">")] = s[1]
    node2parent = {}
    node2order = {}
    with open(os.path.join(taxa_folder, "nodes.dmp")) as nodesfile:
        for line in nodesfile:
            s = map(strip, line.split("|"))
            node2parent[s[0]] = s[1]
            node2order[s[0]] = s[2]
    node2name = {}
    with open(os.path.join(taxa_folder, "names.dmp")) as nodesfile:
        for line in nodesfile:
            s = map(strip, line.split("|"))
            if s[3] == "scientific name":
                node2name[s[0]] = s[1]
    acc2path = {}
    for acc, node in acc2node.items():
        acc2path.setdefault(acc, [])
        parent = node2parent[node]
        while node != parent:
            acc2path[acc].append(node)
            node = parent
            parent = node2parent[node]
    return acc2path, node2order, node2name


def define_common_taxonomy(cluster_dict, node2order, acc2path, node2name, id_perc):
    with open("tax_%i_ref_ITSoneDB.fa" % int(id_perc * 100), "w") as tmp:
        accepted_rank = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        for rep, seq_list in cluster_dict.items():
            if len(seq_list) == 1:
                common_path = acc2path[rep]
                common_path.reverse()
                basic_list = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]
                for node in common_path:
                    if node2order[node] in accepted_rank:
                        index = accepted_rank.index(node2order[node])
                        basic_list[index] = node2order[node][0] + "__" + node2name[node]
                print rep, "singleton", ";".join([i for i in basic_list])
                tmp.write("%s\t%s\n" % (rep, ";".join([i for i in basic_list])))
            else:
                rep_path = acc2path[rep]
                rep_path.reverse()
                path_list = []
                for acc in seq_list:
                    path = acc2path[acc]
                    path_list.reverse()
                    path_list.append(path)
                common_path = []
                for node in rep_path:
                    counts = [path.count(node) for path in path_list]
                    if counts.count(0) == 0:
                        common_path.append(node)
                    else:
                        break
                basic_list = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]
                for node in common_path:
                    if node2order[node] in accepted_rank:
                        index = accepted_rank.index(node2order[node])
                        basic_list[index] = node2order[node][0] + "__" + node2name[node]
                print rep, len(seq_list), ";".join([i for i in basic_list])
                tmp.write("%s\t%s\n" % (rep, ";".join([i for i in basic_list])))


if __name__ == "__main__":
    info_options = options()
    fasta_its1, ncbi_taxa_folder, perc_id = info_options.itsonedb_total_fasta, info_options.taxa_folder, info_options.identity_clustering
    otu_dict, centroidi = clustering_execution(fasta_its1, perc_id)
    path_dict, order_dict, name_dict = generate_taxonomic_info(ncbi_taxa_folder, fasta_its1)
    define_common_taxonomy(otu_dict, order_dict, path_dict, name_dict, perc_id)
