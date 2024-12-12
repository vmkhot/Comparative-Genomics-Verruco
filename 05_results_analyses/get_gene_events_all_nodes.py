## This script will get gene family events at all nodes (including the leaves) and output both with and without the cutoff

import os
import sys
import glob
import pandas as pd
from pathlib import Path

def get_copies_at_node(rec_file):

    """ Function to get the copies of the gene family at a specific node """
    
    rec_file_text = open(rec_file, "r")
    node_info = {}
    for line in rec_file_text:
        line = line.strip()
        if line.startswith("S_"):
            split_table = line.split("\t")
            node = split_table[1]
            dupes_at_node = float(split_table[2])
            trans_at_node = float(split_table[3])
            loss_at_node = float(split_table[4])
            origin_at_node = float(split_table[5])
            copies_at_node = float(split_table[6])
            node_info[node] = [os.path.basename(rec_file),dupes_at_node,trans_at_node,loss_at_node,origin_at_node,copies_at_node]
        # print(node_info)
    return node_info 

def run_node_reconstruction(uml_dir,cutoff):

    """ Function to find all the gene families at each node in the tree """

    uml_dir = Path(uml_dir)
    uml_files = list(uml_dir.glob("*.uml_rec"))
    df_final = pd.DataFrame({'Gene Family':[0], 'Duplications':[0], 'Transfers':[0], 'Loss':[0], 'Originations':[0], 'Copies':[0]})
    # print(df_final)

    for rec_file in uml_files:
        node_info = get_copies_at_node(rec_file)
        df = pd.DataFrame.from_dict(node_info, orient='index')
        df.rename(columns = {0: 'Gene Family', 1: 'Duplications', 2: 'Transfers', 3: 'Loss', 4: 'Originations', 5: 'Copies'}, inplace=True)
        df.reset_index(inplace=True)
        df_final = pd.concat([df_final, df])
        # print(df_final)

    df_final.rename(columns = {'index': 'Node'}, inplace=True)
    print(df_final)
    df_final.to_csv("Events_at_all_nodes.csv",index=False)    # this is a kind of master file with "which genes at which nodes" information to query

    df_node_cut = df_final.loc[(df_final["Duplications"] >= float(cutoff)) | (df_final["Transfers"] >= float(cutoff)) | (df_final["Loss"] >= float(cutoff)) | (df_final["Originations"] >= float(cutoff)) | (df_final["Copies"] >= float(cutoff))]
    df_node_cut.to_csv(f"Events_at_all_nodes_above_{cutoff}.csv",index=False)

    return df_final

def main():
    df = run_node_reconstruction(uml_dir, cutoff)

# %%
if __name__ == "__main__":
    """ To get all nodes: Check the internal node orders. Remember the first and last internal node orders. 
    If the number of internal nodes (branch nodes) are from 16-30, 16 and 30 should be the first and last internal node orders.

    uml_dir : path/to/dir/with/reconciliations
    cutoff : 0-1, the event should be present in this many sampled reconciliations  """

    uml_dir = sys.argv[1]
    cutoff = sys.argv[2]
    main()