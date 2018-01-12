import gseapy as gp
import pandas as pd

from constants import *


def prerank_gsea(df_ranked_genes,
                 report_save_dir,
                 msigdb_file,
                 top_k=100,
                 report_format="png",
                 processes=4,
                 permutation_num=100):
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
    df = df_ranked_genes.sort(['weight'], ascending=[0])

    # Get top k genes
    if top_k > len(df.index):
        top_k = len(df.index)
    df = df[:top_k]

    # Remove header
    df_top_k_wo_header = pd.DataFrame(df.values)

    # Run prerank
    return gp.prerank(rnk=df_top_k_wo_header,
                      gene_sets=msigdb_file,
                      processes=processes,
                      permutation_num=permutation_num,
                      outdir=report_save_dir, format=report_format)