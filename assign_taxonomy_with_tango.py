import os
from string import strip
import argparse
import subprocess
import shlex


def options():
    parser = argparse.ArgumentParser(description="Format ITSoneDB fasta sequences and taxonomy for QIIME. Assumes VSEARCH in in your path.",
                                     prefix_chars="-")
    parser.add_argument("-f", "--itsonedb_total_fasta", type=str,
                        help="file containing all the ITSoneDB fasta sequences",
                        action="store", required=True)
    parser.add_argument("-t", "--taxa_folder", type=str,
                        help="folder containing the NCBI taxonomy data",
                        action="store", required=True)
    parser.add_argument("-p", "--identity_clustering", type=int,
                        help="Identity percent for clustering [default = 0.99",
                        action="store", required=False, default=0.99)
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
    return acc2node, node2order, node2name, node2parent


# NCBI taxonomy pre-processing
# perl preprocess.pl --taxonomy-type NCBI --taxonomy ../../../data_transfer_rel_138/ncbi_taxonomy_data/nodes.dmp  --output ITSoneDB_r138_ref

def tango_excution(cluster_dict, acc2node):
    wd = os.getcwd()
    with open("data_for_tango.csv", "w") as tmp:
        for acc, lista in cluster_dict.items():
            lista = [acc2node[i] for i in lista]
            tmp.write("%s %s\n" % (acc, " ".join(lista)))
    os.chdir("New_TANGO_perl_version")
    cmd = shlex.split("perl tango.pl --taxonomy ITSoneDB_r138_ref --matches %s --output %s" % (
    os.path.join(wd, "data_for_tango.csv"), os.path.join(wd, "tango_recalssification")))
    p = subprocess.Popen(cmd)
    p.wait()
    acc2ass = {}
    with open(os.path.join(wd, "tango_recalssification")) as a:
        for line in a:
            s = map(strip, line.split("\t"))
            acc2ass[s[0]] = s[-1].split(";")[0]
    os.chdir(wd)
    return acc2ass


def define_common_taxonomy(acc2ass, node2order, node2name, id_perc, node2parent):
    with open("tax_%i_ref_ITSoneDB_tango_based.fa" % int(id_perc * 100), "w") as tmp:
        accepted_rank = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        for acc, node in acc2ass.items():
            path = []
            parent = node2parent[node]
            while node != parent:
                path.append(node)
                node = parent
                parent = node2parent[node]
            path.reverse()
            basic_list = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]
            for node in path:
                if node2order[node] in accepted_rank:
                    index = accepted_rank.index(node2order[node])
                    basic_list[index] = node2order[node][0] + "__" + node2name[node]
            tmp.write("%s\t%s\n" % (acc, ";".join([i for i in basic_list])))


if __name__ == "__main__":
    info_options = options()
    fasta_its1, ncbi_taxa_folder, perc_id = info_options.itsonedb_total_fasta, info_options.taxa_folder, info_options.identity_clustering
    otu_dict, centroidi = clustering_execution(fasta_its1, perc_id)
    acc2node_dict, order_dict, name_dict, parent_dict = generate_taxonomic_info(ncbi_taxa_folder, fasta_its1)
    define_common_taxonomy(acc2node_dict, order_dict, name_dict, perc_id, parent_dict)
