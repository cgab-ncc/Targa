import gseapy as gp
import pandas as pd
import numpy as np
import os
from constants import *


def rescale(val, in_min, in_max, out_min, out_max):
    return out_min + (val - in_min) * ((out_max - out_min) / (in_max - in_min))


def prerank_gsea(df_ranked_genes,
                 report_save_dir,
                 msigdb_file,
                 top_k=100,
                 report_format="svg",
                 processes=4,
                 graph_num=100,
                 permutation_num=100,
                 min_size=15,
                 max_size=500):
    """
    This function runs the GSEA prerank function

    REQUIRES:   df_ranked_genes = pandas DataFrame with two columns with headers (1st column = gene, 2nd column = weight)
                top_k = top k number of genes from the df_ranked_genes to be used for enrichment (specify -1 to use all)
                report_save_dir = directory path to where the GSEA report will be saved
                msigdb_file = file path to the msigdb file
                processes = Number of permutations for significance computation
                permutation_num =
    MODIFIES:   nothing
    EFFECTS:    runs the GSEA prerank function and produces reports
    """
    # Sort the dataframe
    df = df_ranked_genes.sort_values(['weight'], ascending=[False])

    # Normalize weight values from -1 to 1
    genes = df.loc[:,'gene'].values.tolist()
    weight_values = df.loc[:, 'weight'].values.tolist()
    in_min = np.min(np.array(weight_values))
    in_max = np.max(np.array(weight_values))
    normalized_weight_values = []
    for i in weight_values:
        normalized_weight_values.append(rescale(val=float(i), in_min=float(in_min), in_max=float(in_max),
                                                out_min=float(-1), out_max=float(1)))
    df = pd.DataFrame({'gene':genes,'weight':normalized_weight_values})

    # Get top k genes
    if top_k != -1:
        if top_k > len(df.index):
            top_k = len(df.index)
        df = df[:top_k]

    # Create the report save dir
    if not os.path.exists(report_save_dir):
        os.makedirs(report_save_dir)

    # Create a rnk file
    with open(report_save_dir + "/" + "temp.rnk", 'w') as file:
        genes = df.loc[:, 'gene'].values.tolist()
        weights = df.loc[:, 'weight'].values.tolist()
        for i in range(0, len(genes)):
            file.write(str(genes[i]) + "\t" + str(weights[i]) + "\n")

    # Load the rnk file
    rank = pd.read_table(report_save_dir + "/" + "temp.rnk", header=None)

    # Run prerank
    return gp.prerank(rnk=rank,
                      gene_sets=msigdb_file,
                      processes=processes,
                      permutation_num=permutation_num,
                      weighted_score_type=0,
                      graph_num=graph_num,
                      min_size=min_size,
                      max_size=max_size,
                      outdir=report_save_dir, format=report_format)