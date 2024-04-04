import numpy as np
import matplotlib.pyplot as plt
import mrcfile
import MDAnalysis as mda
import tools
from tools import *
import os


path = "/media/didymos/Projects/ribosome_tunnels/methionine/"
os.chdir(path)

######################
# Bacteria ribosomes #
######################
bac_id = ['4ybb']

bac_name = ['E.coli']

bac_dict = dict()
bac_dict_long = dict()
i=0
for item in bac_id:
    bac_dict[item] = dict()
    bac_dict_long[item] = dict()
    bac_dict[item]['name'] = bac_name[i]
    bac_dict_long[item]['name'] = bac_name[i]
    i+=1
