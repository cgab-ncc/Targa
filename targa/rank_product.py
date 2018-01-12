from operator import itemgetter
import math


def rank_product(ranked_list_a, ranked_list_b):
    """
    This function calculates the rank products given two tuples

    REQUIRES:	ranked_list_a = [] where each element is the id
                ranked_list_b = [] where each element is the id
                note that both lists ranked_list_a and ranked_list_b are ranked
                (orders are determined as they are passed as this function arguments)
    MODIFIES:	nothing
    EFFECTS:	returns a list of tuples where (id, rank_product)
    """
    if not isinstance(ranked_list_a, list):
        print("ranked_list_a must be a list")
        return
    if not isinstance(ranked_list_b, list):
        print("ranked_list_b must be a list")
        return
    if set(ranked_list_a) != set(ranked_list_b):
        print("Contents of ranked_list_a must and ranked_list_b must be identical")
        return
    rank_products = {} # key = id, value = rank product
    curr_a_rank = 1
    for a in ranked_list_a:
        curr_b_rank = 1
        for b in ranked_list_b:
            if a == b:
                rank_products[a] = math.sqrt(curr_a_rank * curr_b_rank)
                break
            curr_b_rank = curr_b_rank + 1
        curr_a_rank = curr_a_rank + 1
    return sorted(rank_products.items(), key=itemgetter(1), reverse=False)