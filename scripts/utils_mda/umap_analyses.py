from pathlib import Path
from tqdm.auto import tqdm

import numpy as np
import pandas as pd

def get_longest_interactions_indexes(arr):
    nums = [x/2 for x in arr]
    if not nums:
        return []

    nums_set = set(nums)
    longest_subarray = []
    current_subarray = []
    start_index = 0
    end_index = 0

    for i, num in enumerate(nums):
        if num - 1 not in nums_set:
            current_subarray = [num]
            start_index = i
            current_num = num + 1

            while current_num in nums_set:
                current_subarray.append(current_num)
                current_num += 1

            end_index = i + len(current_subarray) - 1

            if end_index - start_index > len(longest_subarray):
                longest_subarray = current_subarray

    return start_index, end_index

def get_longest_interacting_agg_pairs(pep_type, aggreg_times, cutoff):
    cluster_comb_aggreg_longest_times = {
        cluster_pair: {"peptide":[], "timeframe":[]} \
            for cluster_pair in aggreg_times[pep_type].keys()}

    for cluster_pair, d in aggreg_times[pep_type].items():
        if len(d["timeframe"]) > 1:
            longest_indexes = get_longest_interactions_indexes(d["timeframe"])
            x = d["timeframe"]
            cluster_comb_aggreg_longest_times[cluster_pair]["timeframe"] = \
                x[longest_indexes[0]: longest_indexes[1]]
            cluster_comb_aggreg_longest_times[cluster_pair]["peptide"] = \
                d["peptide"][longest_indexes[0]: longest_indexes[1]]
    cluster_comb_aggreg_longest_times = {
        k:v for k, v in cluster_comb_aggreg_longest_times.items() if len(v["timeframe"])>=cutoff
        }
    return cluster_comb_aggreg_longest_times