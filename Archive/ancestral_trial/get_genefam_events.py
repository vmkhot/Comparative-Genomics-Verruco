#! usr/bin/python3

# %%
import pandas as pd
import argparse
# %%
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-f','--infile', nargs='+',help="Not optional path to Copies_at_each_node.csv file")
    parser.add_argument('-e', '--event', nargs='+', help="Not optional. Choose 1 from Duplications,Transfers,Losses,Originations,Copies")
    parser.add_argument('-n', '--node', help="(int) for node, default=ALL")
    parser.add_argument('-g', '--genefam',help="(int) for genefam, default=ALL")
    parser.add_argument('-c', '--cutoff', help="value for boostrap cutoff, default=0.5, this means the gene family has to appear at that node >50% of the bootstraps" )

    parser.print_usage()
    return parser.parse_args()

# %%
def filter_master(master,event,node=None,genefam=None,cutoff=0.5):
    # print specified event for specified gene family at all nodes
    df_master=pd.read_csv(master, header=0, usecols=['Node','Gene_family','Duplications','Transfers','Losses','Originations','Copies'])

    if node is None and genefam is not None:
        df_filter=df_master.loc[(df_master['Gene_family']==f'{genefam}_MSA.fa.ufboot.ale.uml_rec') & (df_master[event] >= float(cutoff)) , ['Node','Gene_family',str(event)]]
        #print(df_master)

    # print all gene families at specified node for specified event
    elif node is not None and genefam is None:
        df_filter=df_master.loc[(df_master['Node']==int(node)) & (df_master[event] >= float(cutoff)) , ['Node','Gene_family',str(event)]]
        #print(df_master)

    #print specified events for all nodes for all gene families
    elif node is None and genefam is None:
        df_filter=df_master.loc[(df_master[event] >= float(cutoff)) , ['Node','Gene_family',str(event)]]
        #print(df_master)

    return df_filter

# %%
def main():
    print('start of script')
    args=parse_args()
    event_in=args.event[0]
    node_in=args.node
    genefam_in=args.genefam
    cutoff_in=args.cutoff
    master=args.infile[0]

    print("file:", master, '\n',
            "out file: Copies_at_each_node_filtered.csv"
            "event:", event_in,'\n',
            "node:", node_in, '\n',
            "gene family:", genefam_in, '\n',
            "cutoff:", cutoff_in,)
    #print(master)
    #df_master=pd.read_csv(master, header=0, usecols=['Node','Gene_family','Duplications','Transfers','Losses','Originations','Copies'])
    #print(df_master)

    df_filter=filter_master(master=master,event=event_in,node=node_in,genefam=genefam_in,cutoff=cutoff_in)
    df_filter.to_csv('./Copies_at_each_node_filtered.csv', index=False , header=True)
    return

main()


