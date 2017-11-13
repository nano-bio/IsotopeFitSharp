using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using DotNetDoctor.csmatio.io;
using DotNetDoctor.csmatio.types;

namespace IsotopeFit
{
    [Obsolete]
    public static class IFDFile
    {
        /// <summary>
        /// Opens an IFD file for further reading.
        /// </summary>
        /// <param name="path">Path to the IFD file.</param>
        /// <returns>Root element of the IFD data structure.</returns>
        public static MLStructure Open(string path)
        {
            MatFileReader reader = new MatFileReader(path);
            return reader.Content["data"] as MLStructure;
        }

        internal static IFData.Spectrum ReadRawData(MLStructure root)
        {
            return new IFData.Spectrum((root["raw_peakdata"] as MLDouble).GetArray());
        }

        internal static ulong ReadStartIndex(MLStructure root)
        {
            return (ulong)(root["startind"] as MLDouble).GetArray()[0][0];
        }

        internal static ulong ReadEndIndex(MLStructure root)
        {
            return (ulong)(root["endind"] as MLDouble).GetArray()[0][0];
        }

        internal static List<IFData.Molecule> ReadMolecules(MLStructure root)
        {
            //root["molecules"]
            //str = (Arr as MLStructure)["peakdata", i] as MLArray;
            //m.PeakData = (str as MLDouble).GetArray();

            MLArray Mol = root["molecules"] as MLArray;

            int num = Mol.N;
            List<IFData.Molecule> Molecules = new List<IFData.Molecule>(num);

            for (int i = 0; i < num; i++)
            {
                Molecules.Add(ReadMolecule(Mol, i));
            }

            return Molecules;
        }

        private static IFData.Molecule ReadMolecule(MLArray molecule, int index)
        {
            IFData.Molecule M = new IFData.Molecule();

            MLArray dataField;

            dataField = (molecule as MLStructure)["peakdata", index] as MLArray;
            M.PeakData = new IFData.Molecule.IsotopeData((dataField as MLDouble).GetArray());

            dataField = (molecule as MLStructure)["name", index] as MLArray;
            M.Name = (dataField as MLChar).GetString(0);

            dataField = (molecule as MLStructure)["minmass", index] as MLArray;
            M.MinMass = (dataField as MLDouble).GetArray()[0][0];
            
            dataField = (molecule as MLStructure)["maxmass", index] as MLArray;
            M.MaxMass = (dataField as MLDouble).GetArray()[0][0];

            dataField = (molecule as MLStructure)["com", index] as MLArray;
            M.CentreOfMass = (dataField as MLDouble).GetArray()[0][0];

            dataField = (molecule as MLStructure)["minind", index] as MLArray;
            M.MinIndex = (ulong)(dataField as MLDouble).GetArray()[0][0];

            dataField = (molecule as MLStructure)["maxind", index] as MLArray;
            M.MaxIndex = (ulong)(dataField as MLDouble).GetArray()[0][0];

            dataField = (molecule as MLStructure)["area", index] as MLArray;
            M.Area = (dataField as MLDouble).GetArray()[0][0];

            dataField = (molecule as MLStructure)["areaerror", index] as MLArray;
            M.AreaError = (dataField as MLSparse).GetArray()[0][0];

            dataField = (molecule as MLStructure)["rootindex", index] as MLArray;
            M.RootIndex = (ulong)(dataField as MLDouble).GetArray()[0][0];

            return M;
        }

        internal static IFData.Calibration ReadCalibration(MLStructure root)
        {
            IFData.Calibration C = new IFData.Calibration();

            MLArray Cal = root["calibration"] as MLArray;

            MLArray dataField;

            dataField = (Cal as MLStructure)["comlist"] as MLArray;
            C.COMList = IFData.Arr2DToVect((dataField as MLDouble).GetArray());

            dataField = (Cal as MLStructure)["massoffsetlist"] as MLArray;
            C.MassOffsetList = IFData.Arr2DToVect((dataField as MLDouble).GetArray());

            dataField = (Cal as MLStructure)["resolutionlist"] as MLArray;
            C.ResolutionList = IFData.Arr2DToVect((dataField as MLDouble).GetArray());

            dataField = (Cal as MLStructure)["massoffsetmethode"] as MLArray;
            C.MassOffsetMethod = (dataField as MLChar).GetString(0);

            dataField = (Cal as MLStructure)["resolutionmethode"] as MLArray;
            C.ResolutionMethod = (dataField as MLChar).GetString(0);

            dataField = (Cal as MLStructure)["massoffsetparam"] as MLArray;
            C.MassOffsetParam = (ulong)(dataField as MLDouble).GetArray()[0][0];

            dataField = (Cal as MLStructure)["resolutionparam"] as MLArray;
            C.ResolutionParam = (ulong)(dataField as MLDouble).GetArray()[0][0];

            dataField = (Cal as MLStructure)["namelist"];
            for (int i = 0; i < dataField.N; i++)
            {
                C.Namelist.Add(((dataField as MLCell).Cells[i] as MLChar).GetString(0));
            }

            C.Shape = ReadShape((Cal as MLStructure)["shape"]);

            return C;
        }

        private static IFData.Calibration.LineShape ReadShape(MLArray sh)
        {
            IFData.Calibration.LineShape S = new IFData.Calibration.LineShape();

            MLArray dataField;

            dataField = (sh as MLStructure)["form"];
            S.Form = (dataField as MLChar).GetString(0);

            dataField = (sh as MLStructure)["breaks"] as MLArray;
            S.Breaks = IFData.Arr2DToVect((dataField as MLDouble).GetArray());

            dataField = (sh as MLStructure)["coefs"] as MLArray;
            S.Coefs = IFData.Arr2DToMatrix((dataField as MLDouble).GetArray());

            dataField = (sh as MLStructure)["pieces"] as MLArray;
            S.Pieces = (ulong)(dataField as MLDouble).GetArray()[0][0];

            dataField = (sh as MLStructure)["order"] as MLArray;
            S.Order = (ulong)(dataField as MLDouble).GetArray()[0][0];

            dataField = (sh as MLStructure)["dim"] as MLArray;
            S.Dim = (ulong)(dataField as MLDouble).GetArray()[0][0];

            return S;
        }

        internal static IFData.BaselineCorr ReadBackgroundCorr(MLStructure root)
        {
            IFData.BaselineCorr BC = new IFData.BaselineCorr();

            MLArray BgCorr = root["bgcorrectiondata"] as MLArray;

            MLArray dataField;

            dataField = (BgCorr as MLStructure)["startmass"] as MLArray;
            BC.StartMass = (dataField as MLDouble).GetArray()[0][0];

            dataField = (BgCorr as MLStructure)["endmass"] as MLArray;
            BC.EndMass = (dataField as MLDouble).GetArray()[0][0];

            dataField = (BgCorr as MLStructure)["ndiv"] as MLArray;
            BC.NDiv = (ulong)(dataField as MLDouble).GetArray()[0][0];

            dataField = (BgCorr as MLStructure)["percent"] as MLArray;
            BC.Percent = (dataField as MLDouble).GetArray()[0][0];

            dataField = (BgCorr as MLStructure)["bgm"] as MLArray;
            BC.XAxis = IFData.Arr2DToVect((dataField as MLDouble).GetArray());

            dataField = (BgCorr as MLStructure)["bgy"] as MLArray;
            BC.YAxis = IFData.Arr2DToVect((dataField as MLDouble).GetArray());

            return BC;
        }
    }
}
