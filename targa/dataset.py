import pandas as pd
from constants import *


def load_data(gene_expression_file, stage_data_file):
    gene_expr_data = GeneExpressionData(gene_expression_file=gene_expression_file)
    stage_data = StageData(stage_data_file=stage_data_file)
    return gene_expr_data, stage_data


class GeneExpressionData:

    def __init__(self, gene_expression_file):
        """
        REQUIRES:   gene_expression_file = string
        MODIFIES:   self.__gene_expression_df
        EFFECTS:    populates the variables in MODIFIES
        """
        if not isinstance(gene_expression_file, str):
            print("gene_expression_file must be a string")
            return
        self.__gene_expression_df = pd.read_table(gene_expression_file, sep="\t", header=0, index_col=0)
        self.__gene_expression_df_t = self.__gene_expression_df.transpose()

    def get_gene_expression_df(self):
        return self.__gene_expression_df

    def get_genes(self):
        return self.__gene_expression_df.iloc[:, 0].index.values.tolist()

    def get_all_sample_ids(self):
        return self.__gene_expression_df.columns.values.tolist()

    def get_tumor_sample_ids(self):
        tumor_samples = []
        for column in self.get_all_sample_ids():
            id = str(column)
            id_components = id.split("-")
            sample_type = str(id_components[-1])
            if sample_type == Dataset.SampleTypes.TUMOR:
                tumor_samples.append(id)
        return tumor_samples

    def get_normal_sample_ids(self):
        normal_samples = []
        for column in self.get_all_sample_ids():
            id = str(column)
            id_components = id.split("-")
            sample_type = str(id_components[-1])
            if sample_type == Dataset.SampleTypes.NORMAL:
                normal_samples.append(id)
        return normal_samples

    def get_metastatic_sample_ids(self):
        metastatic_samples = []
        for column in self.get_all_sample_ids():
            id = str(column)
            id_components = id.split("-")
            sample_type = str(id_components[-1])
            if sample_type == Dataset.SampleTypes.METASTATIC:
                metastatic_samples.append(id)
        return metastatic_samples

    def get_all_patient_ids(self):
        patient_ids = set()
        for sample_id in self.get_all_sample_ids():
            id_components = sample_id.split("-")
            patient_id = "-".join(id_components[:3])
            patient_ids.add(patient_id)
        return list(patient_ids)

    def get_patient_id(self, sample_id):
        id_components = sample_id.split("-")
        patient_id = "-".join(id_components[:3])
        return patient_id

    def get_tumor_id(self, patient_id):
        return patient_id + "-" + Dataset.SampleTypes.TUMOR

    def get_normal_id(self, patient_id):
        return patient_id + "-" + Dataset.SampleTypes.NORMAL

    def get_expression(self, sample_id, gene):
        return self.__gene_expression_df[sample_id][gene]

    def get_expressions(self, sample_ids, gene):
        return self.__gene_expression_df_t.ix[sample_ids][gene]


class StageData:

    def __init__(self, stage_data_file):
        """
        REQUIRES:   stage_data_file = string
        MODIFIES:   self.__stage_data_df
        EFFECTS:    populates the variables in MODIFIES
        """
        if not isinstance(stage_data_file, str):
            print("stage_data_file must be a string")
            return
        self.__stage_data_df = pd.read_csv(stage_data_file)

    def get_df(self):
        return self.__stage_data_df

    def get_sample_ids(self):
        return self.__stage_data_df.loc[:,'sample_id'].values.tolist()

    def get_stage(self, sample_id, type):
        filter_cond = (self.__stage_data_df['sample_id'] == sample_id & self.__stage_data_df['type'] == type)
        return self.__stage_data_df.loc[filter_cond, 'stage'].values.tolist()[0]

    def get_stage_sample_ids(self, stage, type):
        filter_cond = (self.__stage_data_df['stage'] == stage) & (self.__stage_data_df['type'] == type)
        return self.__stage_data_df.loc[filter_cond, 'sample_id'].values.tolist()