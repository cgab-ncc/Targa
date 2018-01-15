import sys
sys.path.append("../")
import targa
import numpy as np

# Files necessary for Targa
gene_expression_file = "/home/jinseoklee/Documents/Projects/Targa/data/tcga/2015/PRAD_gene_expression_2015-02-24"
stage_data_file = "/home/jinseoklee/Documents/Projects/Targa/data/tcga/2015/TCGA_PRAD_Sample_Ids_7th_AJCC_Stage_2015.csv"

# Load data
gene_expression_data, stage_data = targa.load_data(gene_expression_file=gene_expression_file,
                                                   stage_data_file=stage_data_file)

# Build features
df_features = targa.build_late_stage_features(gene_expression_data=gene_expression_data,
                                              stage_data=stage_data,
                                              genes=gene_expression_data.get_genes(),
                                              method=targa.constants.Features.BuildMethods.MEDIAN)

# Rank genes - TOPSIS
criteria_weights = [0.2, 0.2, 0.05, 0.05, 0.2, 0.2, 0.1]
criteria_directions = ['+','+','-','-','+','+','+']
df_topsis_ranked = targa.topsis(df_features=df_features,
                                criteria_weights=criteria_weights,
                                criteria_directions=criteria_directions,
                                verbose=False)

# Save TOPSIS results
df_topsis_ranked.to_csv("TCGA_PRAD_TOPSIS_All_Genes_Similarity_Scores.csv",
                        index=False,
                        columns=['rank', 'gene', 'similarity_score'])

# GSEA - prerank
msigdb_file = "/home/jinseoklee/Documents/Projects/Targa/data/msigdb/msigdb.v6.1.symbols.gmt"
df_ranked_genes = df_topsis_ranked.loc[:,['gene', 'similarity_score']]
df_ranked_genes.columns = ['gene', 'weight'] # rename the 'similarity_score' column to 'weight'
gsea_prerank_results = targa.prerank_gsea(df_ranked_genes=df_ranked_genes,
                                          report_save_dir="/home/jinseoklee/Documents/Projects/Targa/examples/GSEA_TCGA_PRAD_All_MSIGDB",
                                          top_k=-1,
                                          msigdb_file=msigdb_file)
