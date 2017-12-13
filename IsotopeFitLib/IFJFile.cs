using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

using Newtonsoft.Json;
using System.Collections.Specialized;

namespace IsotopeFit
{
    public static class IFJFile
    {
        public static dynamic Open(string path)
        {
            string json;

            using (FileStream f = File.Open(path, FileMode.Open, FileAccess.Read))
            {
                using (StreamReader sr = new StreamReader(f))
                {
                    json = sr.ReadToEnd();
                }
            }

            return JsonConvert.DeserializeObject(json);
        }

        internal static IFData.Spectrum ReadRawData(dynamic rootElement)
        {
            return new IFData.Spectrum(rootElement.raw_x.ToObject<double[]>(), rootElement.raw_y.ToObject<double[]>());
        }

        internal static OrderedDictionary ReadMolecules(dynamic rootElement)
        {
            // read the cluster IDs and remove the mass numbers from the strings and also the [] brackets
            List<string> clusterIDList = rootElement.clusterList.ToObject<List<string>>();

            for (int i = 0; i < clusterIDList.Count; i++)
            {
                clusterIDList[i] = clusterIDList[i].Split(new char[] { '-' })[1].Trim().Replace("[", "").Replace("]", "");

                if (!Char.IsNumber(clusterIDList[i].Last()))
                {
                    clusterIDList[i] += "1";    //TODO: this is stupid and it is not helping much
                }
            }

            OrderedDictionary clusters = new OrderedDictionary(clusterIDList.Count);

            // read data about all clusters
            Dictionary<string, double[][]> clPeakData = rootElement.clusterDist.ToObject<Dictionary<string, double[][]>>();
            Dictionary<string, dynamic> clInfo = rootElement.clusterInfo.ToObject<Dictionary<string, dynamic>>();

            int fails = 0;

            //dynamic test = clInfo["1He1Xe1"];

            for (int i = 0; i < clusterIDList.Count; i++)
            {
                if (clInfo.TryGetValue(clusterIDList[i], out dynamic clInfoEntry) && clPeakData.TryGetValue(clusterIDList[i], out double[][] clPeakDataEntry))
                {
                    IFData.Cluster c = ReadCluster(clInfoEntry, clPeakDataEntry);
                    clusters.Add(clusterIDList[i], c);
                }
                else
                {
                    fails++;
                }
            }

            return clusters;
        }

        private static IFData.Cluster ReadCluster(dynamic clusterInfoElement, double[][] clusterPeakData)
        {
            IFData.Cluster c = new IFData.Cluster
            {
                PeakData = new IFData.Cluster.IsotopeData(clusterPeakData),
                CentreOfMass = clusterInfoElement.centerOfMass.ToObject<double>(),
                MassOffset = clusterInfoElement.massoffset == null ? 0 : clusterInfoElement.massoffset.ToObject<double>(),
                Resolution = clusterInfoElement.resolution == null ? 0 : clusterInfoElement.resolution.ToObject<double>()
            };

            return c;
        }

        internal static IFData.Calibration ReadCalibration(dynamic rootElement)
        {
            IFData.Calibration C = new IFData.Calibration()
            {
                ResolutionMethod = rootElement.calib.resolutionmethode,
                ResolutionParam = rootElement.calib.resolutionparam.ToObject<int>(),
                MassOffsetMethod = rootElement.calib.massoffsetmethode,
                MassOffsetParam = rootElement.calib.massoffsetparam.ToObject<int>(),
                NameList = rootElement.calib.namelist.ToObject<List<string>>(),
                Shape = new IFData.Calibration.LineShape(rootElement.calib.shape.breaks.ToObject<double[]>(), rootElement.calib.shape.coeffs.ToObject<double[][]>())    //TODO: the coeffs maybe have to be reversed
            };

            // we need to cut the numbers out from the name list strings and the [] brackets
            for (int i = 0; i < C.NameList.Count; i++)
            {
                C.NameList[i] = C.NameList[i].Split(new char[] { '-' })[1].Trim().Replace("[", "").Replace("]", "");
            }

            return C;
        }

        internal static IFData.BaselineCorr ReadBackgroundCorr(dynamic rootElement)
        {
            IFData.BaselineCorr bc = new IFData.BaselineCorr()
            {
                XAxis = rootElement.bgcorrection.bgm.ToObject<double[]>(),
                YAxis = rootElement.bgcorrection.bgy.ToObject<double[]>()
            };

            return bc;
        }
    }
}
