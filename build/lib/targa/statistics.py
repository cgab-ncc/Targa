import scipy.stats as stats
import math
from lifelines.statistics import logrank_test


def perform_chi_squared_test(observed_list, expected_list):
    """
    This function performs chi squared test
    REQUIRES:	observed_list = []
                expected_list = []
                observed_list and expected_list are formatted as follows:
                        Outcome1	Outcome2
                Group1	50			5
                Group2	10			20

                Then
                observed_list = [50,5]
                expected_list = [10,20]
    MODIFIES:	nothing
    EFFECTS:	returns test_statistic, p_value
    """
    matrix = [observed_list, expected_list]
    test_statistic, p_value, ddof, expected = stats.chi2_contingency(matrix)
    return test_statistic, p_value


def perform_t_test(data_1, data_2):
    """
    This function performs t-test between datasets data_1 and data_2

    REQUIRES:   data_1 = [] (list)
                data_2 = [] (list)
    MODIFIES:   nothing
    EFFECTS:    performs t-test between the two lists data_1 and data_2
                and returns t_test statistic and p-value
    """
    results = stats.ttest_ind(data_1, data_2)
    return results[0], results[1]


def perform_mann_whitney_u_test(data_1, data_2, alternative):
    """
    This function performs mann whitney u test between datasets data_1 and data_2

    REQUIRES:   data_1 = [] (list)
                data_2 = [] (list)
                alternative = look up scipy.stats.mannwhitneyu doc (https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html)
    MODIFIES:   nothing
    EFFECTS:    performs mann whitney u test between the two lists data_1 and data_2
                and returns test statistic and p-value
    """
    results = stats.mannwhitneyu(data_1, data_2, alternative=alternative)
    return results[0], results[1]


def perform_logrank_test(group_1_status, group_1_duration, group_2_status, group_2_duration, alpha=0.95):
    """
    This function performs logrank test using the lifelines package

    REQUIRES:   group_1_duration = [] where each element is the time/duration
                group_1_status = [] where each element is either 1 ('observed') or 0 ('did not observe')
                group_2_duration = [] where each element is the time/duration
                group_2_status = [] where each element is either 1 ('observed') or 0 ('did not observe')
                alpha = float
    MODIFIES:   nothing
    EFFECTS:    returns the test statistic and p-value
    """
    results = logrank_test(group_1_duration, group_2_duration, group_1_status, group_2_status, alpha=.95)
    return results.test_statistic, results.p_value
