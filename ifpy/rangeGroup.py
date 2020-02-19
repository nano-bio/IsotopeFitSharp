import numpy as np
import clr
import os

clr.AddReference(os.path.abspath('lib/IsotopeFitLib'))
from IsotopeFit import PolyInterpolation, PPInterpolation

# this will create less larger groups or more smaller groups of molecule
# groups belonging together TODO: should be possible to set by GUI
searchrange = 1


class RangeGroup:

    def __init__(
            self,
            minind=None,
            maxind=None,
            minmass=None,
            maxmass=None,
            clusters=None,
            resolution=None,
            resolutionerror=None,
            massoffset=None,
            massoffseterror=None,
            com=None):
        self.minind = minind
        self.maxind = maxind
        self.minmass = minmass
        self.maxmass = maxmass
        self.clusters = clusters
        self.resolution = resolution
        self.resolutionerror = resolutionerror
        self.massoffset = massoffset
        self.massoffseterror = massoffseterror
        self.com = com


def get_calibration_data(x, y, param, method, axis):
    method = str(method).lower()

    if method == "flat":
        np.ones(np.shape(axis)) * param

    elif method == 'polynomial':
        a = PolyInterpolation(x, y, param)
        return a.Evaluate(axis)

    elif method == 'spline':  # called SplineNotAKnot in the C# code, see Michals thesis for more infos
        a = PPInterpolation(x, y, PPInterpolation.PPType.SplineNotAKnot)
        return a.Evaluate(axis)

    elif method == 'pchip':
        a = PPInterpolation(x, y, PPInterpolation.PPType.PCHIP)
        return a.Evaluate(axis)

    else:
        raise Exception("Unknown interpolation method")


def sigma_by_calibration(calibration, mass_axis):
    resolution = resolution_by_calibration(
        calibration,
        mass_axis)
    return mass_axis / resolution * (1 / (2 * np.sqrt(2 * np.log(2))))


def resolution_by_calibration(calibration, mass_axis):
    return get_calibration_data(
        calibration.COMList,
        calibration.ResolutionList,
        calibration.ResolutionParam,
        calibration.ResolutionMethod,
        mass_axis)


def mass_offset_by_calibration(calibration, mass_axis):
    return get_calibration_data(
        calibration.COMList,
        calibration.MassOffsetList,
        calibration.MassOffsetParam,
        calibration.MassOffsetMethod,
        mass_axis)


def find_ranges(clusters, calibration):
    """
    Generation of groups of clusters for calibration fitting
    :param clusters: object of the clusters class
    :param calibration: object of the calibration class
    :return: list of groups of clusters
    :rtype: list
    """

    ranges = []

    import matplotlib.pyplot as plt
    for i in range(clusters.Count):

        # first cluster is added without further evaluation
        if i == 0:
            ranges.append(
                RangeGroup(
                    minmass=clusters[0].PeakData.Mass[0],
                    maxmass=clusters[0].PeakData.Mass[-1],
                    clusters=[clusters[0]],
                    com=clusters[0].CentreOfMass)
            )
            continue

        mass_minus = searchrange * sigma_by_calibration(calibration, clusters[i].CentreOfMass)
        mass_plus = searchrange * sigma_by_calibration(calibration, clusters[i - 1].CentreOfMass)

        if clusters[i].PeakData.Mass[0] - mass_minus <= ranges[-1].maxmass + mass_plus:  # cluster mass range overlaps

            if ranges[-1].maxmass < clusters[i].PeakData.Mass[-1]:
                ranges[-1].maxmass = clusters[i].PeakData.Mass[-1]

            if ranges[-1].minmass > clusters[i].PeakData.Mass[0]:
                ranges[-1].minmass = clusters[i].PeakData.Mass[0]

            ranges[-1].clusters.append(clusters[i])

        else:  # create new mass range
            ranges.append(
                RangeGroup(
                    minmass=clusters[i].PeakData.Mass[0],
                    maxmass=clusters[i].PeakData.Mass[-1],
                    clusters=[clusters[i]]))
    plt.show()
    for range_group in ranges:
        print(len(range_group.clusters))

        # com of group
        arealist = [cluster.Abundance for cluster in range_group.clusters]
        comlist = [cluster.CentreOfMass for cluster in range_group.clusters]
        areasum = sum(arealist)
        comtemp = sum(
            [a * b for a, b in zip(comlist, arealist)]
        )

        if areasum == 0:
            for cluster in range_group.clusters:
                comtemp = + cluster.CentreOfMass
            range_group.com = comtemp / len(range_group.clusters)

        else:
            range_group.com = comtemp / areasum

        # resolution
        range_group.resolution = resolution_by_calibration(calibration, range_group.com)
        range_group.resolutionerror = None

        # mass offset
        range_group.massoffset = mass_offset_by_calibration(calibration, range_group.com)
        range_group.massoffseterror = None

    return ranges
