import os
import sys
import pandas as pd

import matplotlib as mplt
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def GetColorByNormalizedValue(cmap_name, norm_value):
    if norm_value < 0 or norm_value > 1:
        print("ERROR: value " + str(norm_value) + ' does not belong to [0, 1]')
    cmap =  mplt.colormaps[cmap_name] #plt.cm.get_cmap(cmap_name)
    color = cmap(norm_value)
    return mplt.colors.rgb2hex(color[:3])

def ColorByPercentIdentity(pi, config):
    min_pi = config.pi_min
    max_pi = config.pi_max
    fraction = (min(max(pi, min_pi), max_pi) - min_pi) / (max_pi - min_pi)
    if config.cmap_reverse:
        fraction = 1 - fraction
    return GetColorByNormalizedValue(config.cmap, fraction)

def ModifyPos(pos, seq_len, strand):
    if strand == '+':
        return pos
    return seq_len - pos + 1

def PrepareDir(dir_name):
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
