import numpy as np
import pandas as pd

def get_pep_from_res(res, pep_type):
    if "WF1a_only" in pep_type:
        pep_dict = \
            {1: (1, 20),
            2: (21, 40),
            3: (41, 60),
            4: (61, 80),
            5: (81, 100),
            6: (101, 120),
            7: (121, 140),
            8: (141, 160)}
    elif "WF2_only" in pep_type:
        pep_dict = {
            1: (1, 25),
            2: (26, 50),
            3: (51, 75),
            4: (76, 100),
            5: (101, 125),
            6: (126, 150),
            7: (151, 175),
            8: (176, 200)}
    else:
        pep_dict = {1: (1, 20),
                2: (21, 40),
                3: (41, 60),
                4: (61, 80),
                5: (81, 105),
                6: (106, 130),
                7: (131, 155),
                8: (156, 180)}

    val = None
    for k, v in pep_dict.items():
        if int(res) in range(v[0], v[1]+1):
            val = k
    return val