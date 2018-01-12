class Dataset:
    class SampleTypes:
        TUMOR = "01"
        METASTATIC = "06"
        NORMAL = "11"
        ALL = [TUMOR, METASTATIC, NORMAL]

class Features:
    class BuildMethods:
        AVERAGE = "average"
        MEDIAN = "median"
        ALL = [AVERAGE, MEDIAN]