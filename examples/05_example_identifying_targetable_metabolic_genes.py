import sys
sys.path.append("../")
import targa
import numpy as np
import pandas as pd


# Files necessary for Targa
gene_expression_file = "/home/jinseoklee/Documents/Projects/Targa/data/tcga/2015/PRAD_gene_expression_2015-02-24"
stage_data_file = "/home/jinseoklee/Documents/Projects/Targa/data/tcga/2015/TCGA_PRAD_Sample_Ids_7th_AJCC_Stage_2015.csv"
metabolic_genes_file = "/home/jinseoklee/Documents/Projects/Targa/data/genes/metabolic_genes_list.txt"

# Load data
gene_expression_data, stage_data = targa.load_data(gene_expression_file=gene_expression_file,
                                                   stage_data_file=stage_data_file)

# # Build features - for all genes
# df_features = targa.build_late_stage_features(gene_expression_data=gene_expression_data,
#                                               stage_data=stage_data,
#                                               genes=gene_expression_data.get_genes(),
#                                               method=targa.constants.Features.BuildMethods.MEDIAN)
# df_features.to_csv("TCGA_PRAD_All_Genes_7_Features.csv")
#
# # Rank genes - TOPSIS
# criteria_weights = [0.2, 0.2, 0.05, 0.05, 0.2, 0.2, 0.1]
# criteria_directions = ['+','+','-','-','+','+','+']
# df_topsis_ranked = targa.topsis(df_features=df_features,
#                                 criteria_weights=criteria_weights,
#                                 criteria_directions=criteria_directions,
#                                 verbose=False)

df_features = pd.read_csv("TCGA_PRAD_All_Genes_7_Features.csv")
df_topsis_ranked = pd.read_csv("TCGA_PRAD_TOPSIS_All_Genes_Similarity_Scores.csv")

# Get the top 75th percentile genes
similarity_scores = df_topsis_ranked.loc[:,'similarity_score'].values.tolist()
cutoff = np.percentile(similarity_scores, 75)
top_75th_percentile_genes = df_topsis_ranked.loc[df_topsis_ranked['similarity_score'] > cutoff, 'gene'].values.tolist()

# Get features - for the top 75th percentile genes
df_selected_genes_features = df_features.loc[df_features['gene'].isin(top_75th_percentile_genes), :]

# Filter the selected genes by applying the 7 conditions
conditions = (df_selected_genes_features['f1'] > 0) & (df_selected_genes_features['f2'] > 0) & \
             (df_selected_genes_features['f3'] < 0) & (df_selected_genes_features['f4'] < 0) & \
             (df_selected_genes_features['f5'] > 0) & (df_selected_genes_features['f6'] > 0) & \
             (df_selected_genes_features['f7'] > 0)

print("Selected and filtered genes (after applying the 7 conditions)")
df_selected_and_filtered_genes_features = df_selected_genes_features.loc[conditions, :]
selected_and_filtered_genes = df_selected_and_filtered_genes_features.loc[:,'gene'].values.tolist()
print("Number of selected & filtered genes : " + str(len(selected_and_filtered_genes)))

# Load metabolic genes
df_metabolic_genes = pd.read_csv(metabolic_genes_file)
metabolic_genes = [i[0] for i in df_metabolic_genes.values.tolist()]

# Identify metabolic genes in our dataset
common_metabolic_genes = list(set.intersection(set(metabolic_genes), set(gene_expression_data.get_genes())))

# Identify metabolic genes common between top_75th_percentile_genes and common_metabolic_genes
selected_metabolic_genes = list(set.intersection(set(common_metabolic_genes), set(top_75th_percentile_genes)))

# Get features - for the selected and filtered metabolic genes
df_selected_and_filtered_metabolic_genes_features = df_selected_and_filtered_genes_features.loc[df_selected_and_filtered_genes_features['gene'].isin(selected_metabolic_genes), :]
selected_and_filtered_metabolic_genes = df_selected_and_filtered_metabolic_genes_features.loc[:,'gene'].values.tolist()

print("Selected and filtered metabolic genes (after applying the 7 conditions)")
print(selected_and_filtered_metabolic_genes)
print("Number of selected & filtered metabolic genes : " + str(len(selected_and_filtered_metabolic_genes)))

# Get the TOPSIS data of the selected metabolic genes
df_ranked_selected_metabolic_genes = df_topsis_ranked.loc[df_topsis_ranked['gene'].isin(selected_metabolic_genes), :]

# GSEA - prerank
msigdb_file = "/home/jinseoklee/Documents/Projects/Targa/data/msigdb/c2.cgp.v6.1.symbols.gmt"
df_ranked_selected_metabolic_genes = df_ranked_selected_metabolic_genes.loc[:,['gene', 'similarity_score']]
df_ranked_selected_metabolic_genes.columns = ['gene', 'weight'] # rename the 'similarity_score' column to 'weight'
gsea_prerank_results = targa.prerank_gsea(df_ranked_genes=df_ranked_selected_metabolic_genes ,
                                          report_save_dir="/home/jinseoklee/Documents/Projects/Targa/examples/GSEA_TCGA_PRAD_Selected_Metabolic_Genes_Chemical_Genetic_Perturbations",
                                          top_k=-1,
                                          graph_num=150,
                                          msigdb_file=msigdb_file)