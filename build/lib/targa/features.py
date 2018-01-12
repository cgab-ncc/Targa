import pandas as pd
import numpy as np
from constants import *
from dataset import *


def __compute_feature(gene_expression_data, group_1_sample_ids, group_2_sample_ids, gene, method):
    group_1_tumor_sample_expr = np.array(gene_expression_data.get_expressions(sample_ids=group_1_sample_ids, gene=gene))
    group_2_tumor_sample_expr = np.array(gene_expression_data.get_expressions(sample_ids=group_2_sample_ids, gene=gene))
    if (method == Features.BuildMethods.AVERAGE):
        return np.mean(group_1_tumor_sample_expr) - np.mean(group_2_tumor_sample_expr)
    if (method == Features.BuildMethods.MEDIAN):
        return np.median(group_1_tumor_sample_expr) - np.median(group_2_tumor_sample_expr)

def __compute_features(gene_expression_data, stage_data, method, gene):
    # Tumor sample ids
    stage_4_tumor_sample_ids = stage_data.get_stage_sample_ids(stage="4", type="tumor")
    stage_3_tumor_sample_ids = stage_data.get_stage_sample_ids(stage="3", type="tumor")
    stage_2a_tumor_sample_ids = stage_data.get_stage_sample_ids(stage="2A", type="tumor")
    stage_2b_tumor_sample_ids = stage_data.get_stage_sample_ids(stage="2B", type="tumor")
    stage_2_tumor_sample_ids = stage_2a_tumor_sample_ids + stage_2b_tumor_sample_ids

    # Normal sample ids
    stage_4_normal_sample_ids = stage_data.get_stage_sample_ids(stage="4", type="normal")
    stage_3_normal_sample_ids = stage_data.get_stage_sample_ids(stage="3", type="normal")
    stage_2a_normal_sample_ids = stage_data.get_stage_sample_ids(stage="2A", type="normal")
    stage_2b_normal_sample_ids = stage_data.get_stage_sample_ids(stage="2B", type="normal")
    stage_2_normal_sample_ids = stage_2a_normal_sample_ids + stage_2b_normal_sample_ids

    # f1. Stage_4_Tumor_Stat_Summary - Stage_3_Tumor_Stat_Summary
    f1 = __compute_feature(gene_expression_data=gene_expression_data,
                           group_1_sample_ids=stage_4_tumor_sample_ids,
                           group_2_sample_ids=stage_3_tumor_sample_ids,
                           gene=gene, method=method)

    # f2. Stage_3_Tumor_Stat_Summary - Stage_2_Tumor_Stat_Summary
    f2 = __compute_feature(gene_expression_data=gene_expression_data,
                           group_1_sample_ids=stage_3_tumor_sample_ids,
                           group_2_sample_ids=stage_2_tumor_sample_ids,
                           gene=gene, method=method)

    # f3. Stage_4_Normal_Stat_Summary - Stage_3_Normal_Stat_Summary
    f3 = __compute_feature(gene_expression_data=gene_expression_data,
                           group_1_sample_ids=stage_4_normal_sample_ids,
                           group_2_sample_ids=stage_3_normal_sample_ids,
                           gene=gene, method=method)

    # f4. Stage_3_Normal_Stat_Summary - Stage_2_Normal_Stat_Summary
    f4 = __compute_feature(gene_expression_data=gene_expression_data,
                           group_1_sample_ids=stage_3_normal_sample_ids,
                           group_2_sample_ids=stage_2_normal_sample_ids,
                           gene=gene, method=method)

    # f5. Stage_4_Tumor_Stat_Summary - Stage_4_Normal_Stat_Summary
    f5 = __compute_feature(gene_expression_data=gene_expression_data,
                           group_1_sample_ids=stage_4_tumor_sample_ids,
                           group_2_sample_ids=stage_4_normal_sample_ids,
                           gene=gene, method=method)

    # f6. Stage_3_Tumor_Stat_Summary - Stage_3_Normal_Stat_Summary
    f6 = __compute_feature(gene_expression_data=gene_expression_data,
                           group_1_sample_ids=stage_3_tumor_sample_ids,
                           group_2_sample_ids=stage_3_normal_sample_ids,
                           gene=gene, method=method)

    # f7. Stage_2_Tumor_Stat_Summary - Stage_2_Normal_Stat_Summary
    f7 = __compute_feature(gene_expression_data=gene_expression_data,
                           group_1_sample_ids=stage_2_tumor_sample_ids,
                           group_2_sample_ids=stage_2_normal_sample_ids,
                           gene=gene, method=method)

    return f1, f2, f3, f4, f5, f6, f7

def build_late_stage_features(gene_expression_data, stage_data, method=Features.BuildMethods.MEDIAN):
    """
    This function creates a pandas DataFrame of 7 gene expression features focused on late cancer stage

    REQUIRES:   gene_expression_data = an instance of GeneExpressionData
                stage_data = an instance of StageData
                method = string (refer to constants.py, must be one of constants.Features.BuildMethods
    MODIFIES:   nothing
    EFFECTS:    computes 7 features (see below) based on the specified method,
                and returns a pd.DataFrame where there are 7 features (columns) for gene (row)
                f1 = Stage_4_Tumor_Stat_Summary - Stage_3_Tumor_Stat_Summary
                f2 = Stage_3_Tumor_Stat_Summary - Stage_2_Tumor_Stat_Summary
                f3 = Stage_4_Normal_Stat_Summary - Stage_3_Normal_Stat_Summary
                f4 = Stage_3_Normal_Stat_Summary - Stage_2_Normal_Stat_Summary
                f5 = Stage_4_Tumor_Stat_Summary - Stage_4_Normal_Stat_Summary
                f6 = Stage_3_Tumor_Stat_Summary - Stage_3_Normal_Stat_Summary
                f7 = Stage_2_Tumor_Stat_Summary - Stage_2_Normal_Stat_Summary
    """
    if not isinstance(gene_expression_data, GeneExpressionData):
        print("gene_expression_data must be an instance of the class GeneExpressionData")
        return
    if not isinstance(stage_data, StageData):
        print("stage_data must be an instance of the class StageData")
        return
    if not isinstance(method, str):
        print("method must be a string")
        return
    if (method not in set(Features.BuildMethods.ALL)):
        print("method must be one of the constants in Features.BuildMethods")
        return
    data = {'gene','f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7'}
    for gene in gene_expression_data.get_genes():
        f1, f2, f3, f4, f5, f6, f7 = __compute_features(gene_expression_data=gene_expression_data,
                                                        stage_data=stage_data,
                                                        gene=gene)
        data['gene'].append(gene)
        data['f1'].append(f1)
        data['f2'].append(f2)
        data['f3'].append(f3)
        data['f4'].append(f4)
        data['f5'].append(f5)
        data['f6'].append(f6)
        data['f7'].append(f7)

    return pd.DataFrame(data)





