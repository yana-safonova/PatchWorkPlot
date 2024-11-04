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

def ColorByPercentIdentity(cmap, pi, min_pi, max_pi, cmap_reverse):
    fraction = (min(max(pi, min_pi), max_pi) - min_pi) / (max_pi - min_pi)
    if cmap_reverse:
        fraction = 1 - fraction
    return GetColorByNormalizedValue(cmap, fraction)

def rgb2hex(r,g,b):
    return "#{:02x}{:02x}{:02x}".format(r,g,b)

def hex2rgb(hexcode):
    hex_string = hexcode.lstrip('#')
    r = int(hex_string[0:2], 16)
    g = int(hex_string[2:4], 16)
    b = int(hex_string[4:6], 16)
    return (r, g, b)

def ModifyPos(pos, seq_len, strand):
    if strand == '+':
        return pos
    return seq_len - pos + 1

def PrepareDir(dir_name):
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
