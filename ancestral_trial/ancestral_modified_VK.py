# %%
import os
import sys
import glob
import pandas as pd

# %%
def get_copies_at_node(rec_file, node):

    """ Function to get the copies of the gene family at a specific node """

    node = str(node)
    dupes_at_node = 0
    trans_at_node = 0
    loss_at_node = 0
    origin_at_node = 0
    copies_at_node = 0
    rec_file_text = open(rec_file, "r")
    for line in rec_file_text:
        line = line.strip()
        if line.startswith("S_i"):
            split_table = line.split("\t")
            internal_node = split_table[1]
            if internal_node == node:
                dupes_at_node += float(split_table[2])
                trans_at_node += float(split_table[3])
                loss_at_node += float(split_table[4])
                origin_at_node += float(split_table[5])
                copies_at_node += float(split_table[6])

    return [dupes_at_node, trans_at_node, loss_at_node, origin_at_node, copies_at_node]

# %%
def run_node_reconstruction(n1,n2):

    """ Function to find all the gene families at each node in the tree """

    nodes = list(range(n1,n2))
    Gene_families = []
    Nodes = []
    Duplications =[]
    Transfers =[]
    Losses= []
    Originations =[]
    Copies =[]

    for node in nodes:
        for rec_file in glob.glob("*.uml_rec"):
            genefamily_dupe = get_copies_at_node(rec_file,node)[0]
            genefamily_trans = get_copies_at_node(rec_file,node)[1]
            genefamily_loss = get_copies_at_node(rec_file,node)[2]
            genefamily_origin = get_copies_at_node(rec_file,node)[3]
            genefamily_copy = get_copies_at_node(rec_file,node)[4]
            
            Gene_families.append(rec_file)
            Nodes.append(node)
            Duplications.append(genefamily_dupe)
            Transfers.append(genefamily_trans)
            Losses.append(genefamily_loss)
            Originations.append(genefamily_origin)
            Copies.append(genefamily_copy)

    data = {"Node":Nodes, "Gene_family":Gene_families, "Duplications":Duplications, "Transfers":Transfers, "Losses":Losses, "Originations":Originations, "Copies":Copies}
    df = pd.DataFrame.from_dict(data)
    df.to_csv("Copies_at_each_node.csv")    # this is a kind of master file with "which genes at which nodes" information to query

    return df

# %%
def make_node_data(df, cutoff, n1, n2):

    """ Function to find all genes at specific nodes, based on cutoff """

    total_at_node = {}
    nodes = list(range(n1,n2))
    for node in nodes:

        df_node = df[df["Node"]==node]
        total_dupes_at_node = df_node.Duplications.sum()
        total_trans_at_node = df_node.Transfers.sum()
        total_loss_at_node = df_node.Losses.sum()
        total_origin_at_node = df_node.Originations.sum()
        total_copies_at_node = df_node.Copies.sum()
        total_at_node[node] = [total_dupes_at_node , total_trans_at_node , total_loss_at_node , total_origin_at_node , total_copies_at_node]
        #print(total_at_node)
        df_node_cut = df_node.loc[(df_node["Duplications"] >= float(cutoff)) | (df_node["Transfers"] >= float(cutoff)) | (df_node["Losses"] >= float(cutoff)) | (df_node["Originations"] >= float(cutoff)) | (df_node["Copies"] >= float(cutoff))]
        #print(df_node_cut)
        
        if len(df_node_cut) >= 1:
            file_name = ("Node_" + str(node) + "_genes_present.csv")
            df_node_cut.to_csv(file_name)

    df_total = pd.DataFrame.from_dict(total_at_node, orient="index")
    df_total = df_total.reset_index()
    df_total.rename(columns={df_total.columns[0]:"Node",
                             df_total.columns[1]:"Total_Duplications",
                             df_total.columns[2]:"Total_Transfers",
                             df_total.columns[3]:"Total_Losses",
                             df_total.columns[4]:"Total_Originations",
                             df_total.columns[5]:"Total_copies"},
                             inplace = True)

    file_name = "Sum_of_DTLSC_at_each_node.csv"
    df_total.to_csv(file_name)

    return(df_total)

# %%
def clean_up():

    """ Function to clean up directory """

    os.mkdir("Gene_families_at_each_node")
    os.system("mv *present.csv* Gene_families_at_each_node")
    os.mkdir("Total_copies_at_node")
    os.system("mv *.csv Total_copies_at_node")

    print("Gene content at each node in dir. 'Gene_families_at_each_node'")
    print("Total copy at each node in dir. 'Total_copies_at_node'")

# %%
def main():
    df = run_node_reconstruction(n1,n2)
    df_total = make_node_data(df, cutoff, n1, n2)
    clean_up()

# %%
if __name__ == "__main__":
    cutoff = sys.argv[1]
    n1 = int(sys.argv[2])
    n2 = int(sys.argv[3])
    main()


