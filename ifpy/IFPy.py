from System.Collections.Specialized import OrderedDictionary

from IsotopeFit import PyGUI as pg
from IsotopeFit import IFData

def fit(massAxis, signalAxis, peakDataDict):
	# create and fill the .NET dictionary
	dotNetDict = OrderedDictionary()

	for key in peakDataDict.keys():
		cluster = IFData.Cluster();

		peakData = IFData.Cluster.IsotopeData(peakDataDict[key]);
		cluster.PeakData = peakData

		dotNetDict.Add(key, cluster)

	# call the .NET fit function
	pg.TestFit2(dotNetDict)

	# put the results back to a python dictionary
	abundanceDict = {}

	for key in peakDataDict.keys():
		abundanceDict[key] = {}
		abundanceDict[key]["abundance"] = dotNetDict[key].Abundance
		abundanceDict[key]["abundanceError"] = dotNetDict[key].AbundanceError

	return abundanceDict


