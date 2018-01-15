import pandas as pd
import numpy as np
import math
from scipy.spatial import distance
from operator import itemgetter


def topsis(df_features, criteria_weights, criteria_directions, verbose=True):
    """
    REQUIRES:	df_features = pandas DataFrame of 6 columns ('gene', 'f1', 'f2', 'f3', 'f4', 'f5')
                criteria_weights = [] where each element is the weight of criterion in X.
                Note that the i-th element in weight must match up the j-th element in X
                criteria_directions = [] where each element is either '+' or '-'
                '+' indicates that the higher the criterion value, the better
                '-' indicates that the lower the criterion value, the better
    MODIFIES:	nothing
    EFFECTS:	runs the TOPSIS (Tehcnique for Order of Preference by Simiarilty to Ideal Solution)
                and outputs the resulting matrix
                returns [] where each element is a tuple (index, score) - the list is sorted
    """
    # Check the arguments
    if not isinstance(df_features, pd.DataFrame):
        print("df_features must be a pandas DataFrame")
        return
    if not isinstance(criteria_weights, list):
        print("criteria_weights must be a list")
        return
    if abs(sum(criteria_weights) - 1.0) > 0.001:
        print(abs(sum(criteria_weights) - 1.0))
        print("Sum of criteria_weights must be 1.0")
        return
    if not isinstance(criteria_directions, list):
        print("criteria_directions must be a list")
        return
    for d in criteria_directions:
        if d != '+' and d != '-':
            print("criteria_directions values must be either '+' or '-'")
            return
    # Prepare data to run TOPSIS
    genes = df_features.loc[:, 'gene'].values.tolist()
    matrix = df_features.loc[:, df_features.columns != 'gene']
    matrix = matrix.values.tolist()

    # Run TOPSIS
    t = TOPSIS(verbose=verbose)
    results = t.run(matrix=np.matrix(matrix),
                    criteria_weights=criteria_weights,
                    criteria_directions=criteria_directions) # where each element is a tuple (index, score)

    # Prepare results
    data = {'rank':[], 'gene':[], 'similarity_score':[]}
    rank = 1
    for i in results:
        print(i)
        data['rank'].append(rank)
        data['gene'].append(genes[i[0]])
        data['similarity_score'].append(i[1])
        rank = rank + 1
    return pd.DataFrame(data)


class TOPSIS:

    def __init__(self, verbose=True):
        self.__verbose = verbose
        self.matrix = None
        self.pos_ideal_soln = [] # positive ideal solution values
        self.neg_ideal_soln = [] # negative ideal solution values
        self.l2_from_pos_ideal_soln = [] # where each element is the l2 distance from pos_ideal_soln
        self.l2_from_neg_ideal_soln = [] # where each element is the l2 distance from neg_ideal_soln
        self.similarity_scores = {} # where key = index i (from self.matrix) value = similarity score

    def compute_ranking(self):
        """
        This is a helper function for run_topsis

        REQUIRES:   self.similarity_scores is populated
        MODIFIES:   self.alternatives_ranking
        EFFECTS:    ranks the alternatives
                    returns [] where each element is a tuple (index, score)
        """
        sorted_map = sorted(self.similarity_scores.iteritems(), key=itemgetter(1), reverse=True)
        if self.__verbose:
            print("Sorted/Ranked Similarity Scores")
            print("Rank,\tindex,\tsimilarity score")
        rank = 1
        ranking = [] # where each element is a tuple (index, score)
        for s in sorted_map:
            ranking.append((int(s[0])-1,float(s[1])))
            if self.__verbose:
                print(str(rank) + ",\t" + str(int(s[0])-1) + ",\t" + str(s[1]))
            rank = rank + 1
        if self.__verbose:
            print("\n")
        return ranking

    def compute_similarity_scores(self):
        """
        This is a helper function for run_topsis

        REQUIRES:   self.l2_from_pos_ideal_soln is populated
                    self.l2_from_neg_ideal_soln is populated
        MODIFIES:   self.similarity_scores
        EFFECTS:    calculates the similarity scores for each alternative
        """
        if self.__verbose:
            print("Similarity Scores (Unordered")
            print("alternative index i (+1),\tscore")
        for i in range(0, len(self.l2_from_pos_ideal_soln)):
            score = self.l2_from_neg_ideal_soln[i] / (self.l2_from_pos_ideal_soln[i] + self.l2_from_neg_ideal_soln[i])
            self.similarity_scores[str(i + 1)] = score
            if self.__verbose:
                print(str(i + 1) + ",\t" + str(score))
        if self.__verbose:
            print("\n")

    def compute_l2_distances(self):
        """
        This is a helper function for run_topsis

        REQUIRES:   self.matrix is populated
                    self.pos_ideal_soln is populated
                    self.neg_ideal_soln is populated
        MODIFIES:   self.l2_from_pos_ideal_soln
                    self.l2_from_neg_ideal_soln
        EFFECTS:    self.l2_from_pos_ideal_soln and
                    self.l2_from_neg_ideal_soln
                    are populated
        """
        for i in self.matrix.A:
            self.l2_from_pos_ideal_soln.append(distance.euclidean(i, self.pos_ideal_soln))
            self.l2_from_neg_ideal_soln.append(distance.euclidean(i, self.neg_ideal_soln))

    def find_pos_neg_ideal_solutions(self, criteria_directions):
        """
        This is a helper function for run_topsis

        REQUIRES:	self.matrix is populated
        MODIFIES:	self.neg_ideal_solution and self.pos_ideal_solution
        EFFECTS:	identifies negative (worst) and positive (best) idea solutions
        """
        i_criterion = 0
        for i in self.matrix.T.A:
            if criteria_directions[i_criterion] == '+':
                # The higher the better
                self.pos_ideal_soln.append(max(i))
                self.neg_ideal_soln.append(min(i))
            else:
                # The lower the better
                self.pos_ideal_soln.append(min(i))
                self.neg_ideal_soln.append(max(i))
            i_criterion = i_criterion + 1

    def normalize_matrix(self, matrix, criteria_weights):
        """
        This is a helper function for run_topsis

        REQUIRES:	matrix = np.matrix
                    criteria_weight = []
        MODIFIES:	self.matrix
        EFFECTS:	populates self.matrix with the normalized values
        """
        # Obtain the summation of (each value)^2, grouped by criterion
        weights_raw = []
        for i in matrix.T.A:
            curr_criterion_weight = 0
            for j in i:
                curr_criterion_weight = curr_criterion_weight + math.pow(j, 2)
            weights_raw.append(curr_criterion_weight)
        # Take the square root of each weight
        weights_sq_root = []
        for w in weights_raw:
            weights_sq_root.append(math.sqrt(w))
        # Perform normalization, including multiplying by the criterion weight
        X_normalized = []
        i_criterion = 0
        for i in matrix.T.A:
            curr_criterion_normalized = []
            for j in i:
                j_normalized = (j / weights_sq_root[i_criterion]) * criteria_weights[i_criterion]
                curr_criterion_normalized.append(j_normalized)
            X_normalized.append(curr_criterion_normalized)
            i_criterion = i_criterion + 1
        # Transpose X_normalized
        X_normalized = np.matrix(X_normalized)
        X_normalized = X_normalized.T
        # Copy X_normalized to self.matrix
        self.matrix = []
        for i in X_normalized.A:
            self.matrix.append(i)
        self.matrix = np.matrix(self.matrix)

    def run(self, matrix, criteria_weights, criteria_directions):
        """
        REQUIRES:	matrix = np.matrix
                    criteria_weights = [] where each element is the weight of criterion in X.
                    Note that the i-th element in weight must match up the j-th element in X
                    criteria_directions = [] where each element is either '+' or '-'
                    '+' indicates that the higher the criterion value, the better
                    '-' indicates that the lower the criterion value, the better
        MODIFIES:	nothing
        EFFECTS:	runs the TOPSIS (Tehcnique for Order of Preference by Simiarilty to Ideal Solution)
                    and outputs the resulting matrix
                    returns [] where each element is a tuple (index, score) - the list is sorted
        """
        # Check the arguments
        if not isinstance(matrix, np.matrix):
            print("X must be an evaluation matrix")
            return
        if not isinstance(criteria_weights, list):
            print("criteria_weights must be a list")
            return
        if abs(sum(criteria_weights) - 1.0) > 0.001:
            print(abs(sum(criteria_weights) - 1.0))
            print("Sum of criteria_weights must be 1.0")
            return
        if not isinstance(criteria_directions, list):
            print("criteria_directions must be a list")
            return
        for d in criteria_directions:
            if d != '+' and d != '-':
                print("criteria_directions values must be either '+' or '-'")
                return
        if self.__verbose:
            print("## Started running TOPSIS ##")
        # Step 1. Normalize the matrix
        self.normalize_matrix(matrix=matrix, criteria_weights=criteria_weights)
        # Step 2. Compute positive and negative ideal solutions
        self.find_pos_neg_ideal_solutions(criteria_directions=criteria_directions)
        # Step 3. Compute L2 distance from positive and negative ideal solutions
        self.compute_l2_distances()
        # Step 4. Compute similarity scores for each alternative
        self.compute_similarity_scores()
        # Step 5. Rank alternatives
        ranking = self.compute_ranking()
        if self.__verbose:
            print("## Finished running TOPSIS ##\n")
        return ranking

# Test
# X = []
# X.append([0.365, 0.00000001, 0.00005])
# X.append([0.165, 0.00000005, 0.00001])
# X.append([0.335, 0.00000006, 0.00003])
# X.append([0.065, 0.00000021, 0.00009])
# X.append([0.765, 0.00000011, 0.00002])
# X.append([-0.365, 0.00000031, 0.00006])
#
# a = 0.33333333
# criteria_weights = [a,a,a]
# criteria_directions = ['+','-','-']
#
# t = TOPSIS()
# t.run(matrix=np.matrix(X),
# 	  criteria_weights=criteria_weights,
# 	  criteria_directions=criteria_directions)