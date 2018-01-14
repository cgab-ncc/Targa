import sys
sys.path.append("../")
import targa
import numpy as np
import pandas as pd


# Files necessary for Targa
# gene_expression_file = "/home/jinseoklee/Documents/Projects/Targa/data/samples/sample_gene_expression_data"
# stage_data_file = "/home/jinseoklee/Documents/Projects/Targa/data/samples/sample_stage_data.csv"
gene_expression_file = "/home/jinseoklee/Documents/Projects/Targa/data/tcga/PRAD_gene_expression_HiSeqV2_2017-10-13"
stage_data_file = "/home/jinseoklee/Documents/Projects/Targa/data/tcga/TCGA_PRAD_Sample_Ids_7th_AJCC_Stage.csv"
metabolic_genes_file = "/home/jinseoklee/Documents/Projects/Targa/data/genes/metabolic_genes_list.txt"

# Load data
gene_expression_data, stage_data = targa.load_data(gene_expression_file=gene_expression_file,
                                                   stage_data_file=stage_data_file)
df_metabolic_genes = pd.read_csv(metabolic_genes_file)
metabolic_genes = [i[0] for i in df_metabolic_genes.values.tolist()]

# Identify metabolic genes in our dataset
common_metabolic_genes = list(set.intersection(set(metabolic_genes), set(gene_expression_data.get_genes())))

# Build features - for only metabolic genes
df_features = targa.build_late_stage_features(gene_expression_data=gene_expression_data,
                                              stage_data=stage_data,
                                              genes=common_metabolic_genes,
                                              method=targa.constants.Features.BuildMethods.MEDIAN)

# Rank genes - TOPSIS
criteria_weights = [0.25, 0.2, 0.25, 0.2, 0.1]
criteria_directions = ['+','+','+','+','+'] # the higher the feature value, the closer a given gene is to the ideal solution
df_topsis_ranked = targa.topsis(df_features=df_features,
                                criteria_weights=criteria_weights,
                                criteria_directions=criteria_directions,
                                verbose=False)

# Save TOPSIS results
df_topsis_ranked.to_csv("TCGA_PRAD_TOPSIS_Metabolic_Genes_Similarity_Scores.csv", index=False, columns=['rank', 'gene', 'similarity_score'])

# GSEA - prerank
msigdb_file = "/home/jinseoklee/Documents/Projects/Targa/data/msigdb/msigdb.v6.1.symbols.gmt"
df_ranked_genes = df_topsis_ranked.loc[:,['gene', 'similarity_score']]
df_ranked_genes.columns = ['gene', 'weight'] # rename the 'similarity_score' column to 'weight'
gsea_prerank_results = targa.prerank_gsea(df_ranked_genes=df_ranked_genes,
                                          report_save_dir="/home/jinseoklee/Documents/Projects/Targa/examples/GSEA_TCGA_PRAD_Metabolic_Genes_All_MSIGDB_Genes",
                                          top_k=-1,
                                          msigdb_file=msigdb_file)
