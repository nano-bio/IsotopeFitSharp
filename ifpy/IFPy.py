"""
this file will cover the work progress how it should be done in GUI but without GUI.
This will be done by some basic examples:
- import data
- correct offset
- make peak shape
- import molecules
- create calibration
- get fit data
"""
# python modules
import os
import clr
# from understanding_clusters import understanding_clusters
import rangeGroup
import matplotlib.pyplot as plt

# C# modules
try:
    clr.AddReference(os.path.abspath('lib/IsotopeFitLib'))
except Exception as e:
    print("Did you use python 64bit?")
from IsotopeFit import Workspace

# create workspace
w = Workspace('testfile.ifd')

# understanding_clusters(w.Clusters)
range_groups = rangeGroup.find_ranges(w.Clusters, w.Calibration)
print(len(range_groups))

for range_group in range_groups:
    mass
    plt.plot()
plt.show()

# load h5 file and save mass and signal data to workspace
# w = load_mass_signal_data(w)

# generate baseline by the two input numbers and correct baseline
# w = generate_and_correct_baseline(w)

# create cluster list


# save cluster list to workspace

print("\n----------------------\nEnd of code in IFPy.py\n----------------------\n")
