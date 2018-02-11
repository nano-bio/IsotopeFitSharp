# from System.Collections.Generic import Dictionary
# from System import String

# params = {"one": "1","two":"2"}   

# dict1 = Dictionary[String, String]()
# for k, v in param.iterkeys():
#    dict1[k] = v

# ret = Hello_Dict(dict1)


import sys
import os
import clr

from System.Collections.Specialized import OrderedDictionary
from System import String, Int32

from IsotopeFit import PyGUI as pg
from IsotopeFit import IFData

def Fit(d):
	dic = OrderedDictionary()

	for k in d.keys():
		c = IFData.Cluster();
		c.Name = k
		c.CentreOfMass = d[k]["mass"]

		isot = IFData.Cluster.IsotopeData(d[k]["isot"]);
		c.PeakData = isot

		dic[k] = c

	pg.TestFit2(dic)

	#return dic


