![alt text](Targa_logo.png)
[![GitHub version](https://badge.fury.io/gh/cgab-ncc%2FTarga.svg)](http://badge.fury.io/gh/cgab-ncc%2FTarga)
[![python 2.7](https://img.shields.io/badge/python-2.7-blue.svg)](https://img.shields.io/badge/python-3.6-blue.svg)
[![Open Source Love](https://badges.frapsoft.com/os/mit/mit.svg?v=102)](https://github.com/ellerbrock/open-source-badge/)

## Targa (TARgetable Genes Analysis)
Targa is a python package for identifying targetable genes in cancer using gene expression data.<br>

## Dependencies
Targa requires:
* Python 2.7
* gseapy (>= 0.9.3)
* numpy
* scipy
* pandas

## Installation
To install Targa, run the install command:
```
pip install targa
```

Alternatively, you can download the Github repository and run the install command:
```
python setup.py install
```

## Usage
Here are example python codes for using Targa<br>

#### 1. Load data
```
import targa as tg
import pandas as pd

gene_expression_file = "sample_gene_expression_file.csv"
clinical_data_file = "sample_clinical_data_file.csv"

gene_expression_data, stage_data = tg.load_data(gene_expression_file=gene_expression_file,
                                                stage_data_file=stage_data_file)
```

#### 2. Build features
```
df_features = targa.build_late_stage_features(gene_expression_data=gene_expression_data,
                                              stage_data=stage_data,
                                              genes=gene_expression_data.get_genes(),
                                              method=tg.constants.Features.BuildMethods.MEDIAN)
```

#### 3. Rank genes by TOPSIS
```
# The higher the feature value, the closer a given gene is to the ideal solution
criteria_directions = ['+','+','+','+','+']

# Late stage focused
criteria_weights = [0.25, 0.2, 0.25, 0.2, 0.1]

df_topsis_ranked = targa.topsis(df_features=df_features,
                                criteria_weights=criteria_weights,
                                criteria_directions=criteria_directions,
                                verbose=False)
```

#### 4. Run GSEA prerank()
```
msigdb_file = "msigdb.v6.1.symbols.gmt"
df_ranked_genes = df_topsis_ranked.loc[:,['gene', 'similarity_score']]

# Rename the 'similarity_score' column to 'weight'
df_ranked_genes.columns = ['gene', 'weight'] 
gsea_prerank_results = targa.prerank_gsea(df_ranked_genes=df_ranked_genes,
                                          report_save_dir="gsea_prerank",
                                          top_k=-1,
                                          msigdb_file=msigdb_file)
```