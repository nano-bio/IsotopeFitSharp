searchrange = 1


class RangeGroup:

    def __init__(
            self,
            minind=None,
            maxind=None,
            minmass=None,
            maxmass=None,
            molecules=None,
            resolution=None,
            resolutionerror=None,
            massoffset=None,
            massoffseterror=None,
            com=None):
        self.minind = minind
        self.maxind = maxind
        self.minmass = minmass
        self.maxmass = maxmass
        self.molecules = molecules
        self.resolution = resolution
        self.resolutionerror = resolutionerror
        self.massoffset = massoffset
        self.massoffseterror = massoffseterror
        self.com = com


def find_ranges(molecules, calibration):
    """
    Generation of groups of molecules for calibration fitting
    :param molecules: object of the Molecules class
    :param calibration: object of the calibration class
    :return: list of groups of molecules
    :rtype: list
    """

    ranges = []
    molecules = molecules.mol_matrix
    for i in range(len(molecules)):

        # first molecule is added without further evaluation
        if i == 0:
            ranges.append(
                RangeGroup(
                    minind=molecules[0]['minind'],
                    maxind=molecules[0]['maxind'],
                    minmass=molecules[0]['minmass'],
                    maxmass=molecules[0]['maxmass'],
                    molecules=[molecules[0]]))
            continue

        mass_minus = searchrange * sigma_by_calibration(calibration, molecules[i]['com'])
        mass_plus = searchrange * sigma_by_calibration(calibration, molecules[i - 1]['com'])

        if molecules[i]['minmass'] - mass_minus <= ranges[-1].maxmass + mass_plus:  # molecule mass range overlaps

            if ranges[-1].maxind < molecules[i]['maxind']:
                ranges[-1].maxind = molecules[i]['maxind']
                ranges[-1].maxmass = molecules[i]['maxmass']

            if ranges[-1].minind > molecules[i]['minind']:
                ranges[-1].minind = molecules[i]['minind']
                ranges[-1].minmass = molecules[i]['minmass']

            ranges[-1].molecules.append(molecules[i])

        else:  # create new massrange
            ranges.append(
                RangeGroup(
                    minind=molecules[i]['minind'],
                    maxind=molecules[i]['maxind'],
                    minmass=molecules[i]['minmass'],
                    maxmass=molecules[i]['maxmass'],
                    molecules=[molecules[i]]))

    for range_group in ranges:

        # com of group
        arealist = [molecule['area'] for molecule in range_group.molecules]
        comlist = [molecule['com'] for molecule in range_group.molecules]
        areasum = sum(arealist)
        comtemp = sum(
            [a * b for a, b in zip(comlist, arealist)]
        )

        if areasum == 0:
            for molecule in range_group.molecules:
                comtemp = + molecule['com']
            range_group.com = comtemp / len(range_group.molecules)

        else:
            range_group.com = comtemp / areasum

        # resolution
        range_group.resolution = resolution_by_calibration(calibration, range_group.com)
        range_group.resolutionerror = None

        # mass offset
        range_group.massoffset = mass_offset_by_calibration(calibration, range_group.com)
        range_group.massoffseterror = None

    return ranges
