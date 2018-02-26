from System.Collections.Specialized import OrderedDictionary

from IsotopeFit import PyGUI as pg
from IsotopeFit import IFData

def fit(massAxis, signalAxis,
		peakDataDict,
		resolutionPPBreaks, resolutionCoeffs,
		peakShapeBreaks, peakShapeCoeffs,
		dmSearchRange, dmFwhmRange):
	"""Function for fitting the spectrum from Python GUI.

	Args:
		massAxis			(list): Corrected mass axis.
		signalAxis			(list): Corrected signal axis.
		peakDataDict		(dict of 2d lists): Structure containing the isotope data for each cluster, with the following structure: { "cluster1":[ [mass1, abundance1], [mass2, abundance2], ...], "cluster2":[...]}
		resolutionPPBreaks	(list): If the resolution fit was a partial polynomial, this contains the PP breaks. Otherwise None.
		resolutionCoeffs	(2d list): Resolution fit coefficients for each partial polynomial. If the root list contains only a single child list, a polynomial interpolation is assumed.
		peakShapeBreaks		(list): Breaks of the line shape fit partial polynomials.
		peakShapeCoeffs		(2d list): Coefficients of the line shape partial polynomials.
		dmSearchRange		(float64): Search range for the design matrix build function.
		dmFwhmRange			(float64): FWHM range for the design matrix build function.

	Returns:
		dict: Dictionary with cluster IDs as keys, containing small 2-key dictionaries with "abundance" and "abundanceError" data for each cluster.
	
	"""
	# create and fill the .NET dictionary
	clusters = OrderedDictionary()

	for key in peakDataDict.keys():
		cluster = IFData.Cluster();

		peakData = IFData.Cluster.IsotopeData(peakDataDict[key]);
		cluster.PeakData = peakData

		clusters.Add(key, cluster)

	# call the .NET fit function
	pg.Fit(massAxis, signalAxis, clusters, resolutionPPBreaks, resolutionCoeffs, peakShapeBreaks, peakShapeCoeffs, dmSearchRange, dmFwhmRange)

	# put the results back to a python dictionary
	abundanceDict = {}

	for key in peakDataDict.keys():
		abundanceDict[key] = {}
		abundanceDict[key]["abundance"] = clusters[key].Abundance
		abundanceDict[key]["abundanceError"] = clusters[key].AbundanceError

	return abundanceDict


