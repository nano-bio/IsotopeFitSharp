using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using DotNetDoctor.csmatio.io;
using DotNetDoctor.csmatio.types;

namespace IsotopeFit
{
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

        internal static int ReadStartIndex(MLStructure root)
        {
            return (int)(root["startind"] as MLDouble).GetArray()[0][0];
        }

        internal static int ReadEndIndex(MLStructure root)
        {
            return (int)(root["endind"] as MLDouble).GetArray()[0][0];
        }

        internal static OrderedDictionary ReadMolecules(MLStructure root)
        {
            MLArray cArr = root["molecules"] as MLArray;

            int num = cArr.N;
            OrderedDictionary clusters = new OrderedDictionary();

            for (int i = 0; i < num; i++)
            {
                IFData.Cluster c = ReadCluster(cArr, i);
                clusters.Add(c.Name, c);
            }

            return clusters;
        }

        private static IFData.Cluster ReadCluster(MLArray molecule, int index)
        {
            IFData.Cluster M = new IFData.Cluster();

            MLArray dataField;

            dataField = (molecule as MLStructure)["peakdata", index] as MLArray;
            M.PeakData = new IFData.Cluster.IsotopeData((dataField as MLDouble).GetArray());

            dataField = (molecule as MLStructure)["name", index] as MLArray;
            M.Name = (dataField as MLChar).GetString(0);

            //dataField = (molecule as MLStructure)["minmass", index] as MLArray;
            //M.MinMass = (dataField as MLDouble).GetArray()[0][0];
            
            //dataField = (molecule as MLStructure)["maxmass", index] as MLArray;
            //M.MaxMass = (dataField as MLDouble).GetArray()[0][0];

            dataField = (molecule as MLStructure)["com", index] as MLArray;
            M.CentreOfMass = (dataField as MLDouble).GetArray()[0][0];

            //dataField = (molecule as MLStructure)["minind", index] as MLArray;
            //M.MinIndex = (int)(dataField as MLDouble).GetArray()[0][0];

            //dataField = (molecule as MLStructure)["maxind", index] as MLArray;
            //M.MaxIndex = (int)(dataField as MLDouble).GetArray()[0][0];

            //dataField = (molecule as MLStructure)["area", index] as MLArray;
            //M.Area = (dataField as MLDouble).GetArray()[0][0];

            //dataField = (molecule as MLStructure)["areaerror", index] as MLArray;
            //M.AreaError = (dataField as MLSparse).GetArray()[0][0];

            //dataField = (molecule as MLStructure)["rootindex", index] as MLArray;
            //M.RootIndex = (int)(dataField as MLDouble).GetArray()[0][0];

            return M;
        }

        internal static IFData.Calibration ReadCalibration(MLStructure root)
        {
            IFData.Calibration C = new IFData.Calibration();

            MLArray Cal = root["calibration"] as MLArray;

            MLArray dataField;

            dataField = (Cal as MLStructure)["comlist"] as MLArray;
            C.COMList = IFData.Arr2DTo1D((dataField as MLDouble).GetArray());

            dataField = (Cal as MLStructure)["massoffsetlist"] as MLArray;
            C.MassOffsetList = IFData.Arr2DTo1D((dataField as MLDouble).GetArray());

            dataField = (Cal as MLStructure)["resolutionlist"] as MLArray;
            C.ResolutionList = IFData.Arr2DTo1D((dataField as MLDouble).GetArray());

            dataField = (Cal as MLStructure)["massoffsetmethode"] as MLArray;
            C.MassOffsetMethod = (dataField as MLChar).GetString(0);

            dataField = (Cal as MLStructure)["resolutionmethode"] as MLArray;
            C.ResolutionMethod = (dataField as MLChar).GetString(0);

            dataField = (Cal as MLStructure)["massoffsetparam"] as MLArray;
            C.MassOffsetParam = (int)(dataField as MLDouble).GetArray()[0][0];

            dataField = (Cal as MLStructure)["resolutionparam"] as MLArray;
            C.ResolutionParam = (int)(dataField as MLDouble).GetArray()[0][0];

            dataField = (Cal as MLStructure)["namelist"];
            for (int i = 0; i < dataField.N; i++)
            {
                C.NameList.Add(((dataField as MLCell).Cells[i] as MLChar).GetString(0));
            }

            C.Shape = ReadShape((Cal as MLStructure)["shape"]);

            return C;
        }

        private static IFData.Calibration.LineShape ReadShape(MLArray sh)
        {
            IFData.Calibration.LineShape S = new IFData.Calibration.LineShape();

            MLArray dataField;

            //dataField = (sh as MLStructure)["form"];
            //S.Form = (dataField as MLChar).GetString(0);

            dataField = (sh as MLStructure)["breaks"] as MLArray;
            S.Breaks = IFData.Arr2DTo1D((dataField as MLDouble).GetArray());

            dataField = (sh as MLStructure)["coefs"] as MLArray;

            double[][] coefs = (dataField as MLDouble).GetArray();

            for (int i = 0; i < coefs.Length; i++)
            {
               coefs[i] = coefs[i].Reverse().ToArray();
            }

            S.Coeffs = IFData.Arr2DToMatrix(coefs);

            //dataField = (sh as MLStructure)["pieces"] as MLArray;
            //S.Pieces = (int)(dataField as MLDouble).GetArray()[0][0];

            //dataField = (sh as MLStructure)["order"] as MLArray;
            //S.Order = (int)(dataField as MLDouble).GetArray()[0][0];

            //dataField = (sh as MLStructure)["dim"] as MLArray;
            //S.Dim = (int)(dataField as MLDouble).GetArray()[0][0];

            return S;
        }

        internal static IFData.BaselineCorr ReadBackgroundCorr(MLStructure root)
        {
            IFData.BaselineCorr BC = new IFData.BaselineCorr();

            MLArray BgCorr = root["bgcorrectiondata"] as MLArray;

            MLArray dataField;

            //dataField = (BgCorr as MLStructure)["startmass"] as MLArray;
            //BC.StartMass = (dataField as MLDouble).GetArray()[0][0];

            //dataField = (BgCorr as MLStructure)["endmass"] as MLArray;
            //BC.EndMass = (dataField as MLDouble).GetArray()[0][0];

            //dataField = (BgCorr as MLStructure)["ndiv"] as MLArray;
            //BC.NDiv = (int)(dataField as MLDouble).GetArray()[0][0];

            //dataField = (BgCorr as MLStructure)["percent"] as MLArray;
            //BC.Percent = (dataField as MLDouble).GetArray()[0][0];

            dataField = (BgCorr as MLStructure)["bgm"] as MLArray;
            BC.XAxis = IFData.Arr2DTo1D((dataField as MLDouble).GetArray());

            dataField = (BgCorr as MLStructure)["bgy"] as MLArray;
            BC.YAxis = IFData.Arr2DTo1D((dataField as MLDouble).GetArray());

            return BC;
        }
    }
}
