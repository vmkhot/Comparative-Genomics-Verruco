# %%
import pandas as pd
import argparse
# %%
def parse_args():
    parser = argparse.ArgumentParser(prog="./get_genefam_events.py")

    parser.add_argument('-e', '--event', nargs='+')
    parser.add_argument('-n', '--node', help="(int) for node, default=ALL")
    parser.add_argument('-g', '--genefam',help="(int) for genefam, default=ALL")
    parser.add_argument('-c', '--cutoff', help="value for boostrap cutoff, default=0.5, this means the gene family has to appear at that node >50% of the time" )
    parser.add_argument('-f','--infile', nargs='+',help="path to Copies_at_each_node.csv file")

    parser.print_usage()
    return parser.parse_args()

# %%
def filter_master(event,node=None,genefam=None,cutoff=0.5):
    # print specified event for specified gene family at all nodes 
    if node is None and genefam is not None:
        df_filter=df_master.loc[(df_master['Gene_family']==f'{genefam}_MSA.fa.ufboot.ale.uml_rec') & (df_master[event] >= cutoff) , ['Node','Gene_family',str(event)]]
        print(df_filter)
    # print all gene families at specified node for specified event
    elif node is not None and genefam is None:
        df_filter=df_master.loc[(df_master['Node']==node) & (df_master[event] >= cutoff) , ['Node','Gene_family',str(event)]]
        print(df_filter)
    #print specified events for all nodes for all gene families
    elif node is None and genefam is None:
        df_filter=df_master.loc[(df_master[event] >= cutoff) , ['Node','Gene_family',str(event)]]
        print(df_filter)
    return df_filter

# %%
def main():
    args=parse_args()
    event_in=args.event
    node_in=int(args.node)
    genefam_in=int(args.genefam)
    cutoff_in=float(args.cutoff)

    master=args.infile
    df_master=pd.read_csv(master, header=0, usecols=['Node','Gene_family','Duplications','Transfers','Losses','Originations','Copies'])
    
    filter_master(event=event_in,node=node_in,genefam=genefam_in,cutoff=cutoff_in)
    return



