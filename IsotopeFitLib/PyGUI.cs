using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using IsotopeFit;

namespace IsotopeFit
{
    public static class PyGUI
    {
        public static void Fit(
            double[] massAxis, double[] signalAxis,    // spectral data
            Dictionary<string, IFData.Cluster> clusters, List<string> clusterIDList,    // cluster data
            double[] resolutionPPBreaks, double[][] resolutionFitCoeffs,    // resolution fit data
            double[] peakShapeBreaks, double[][] peakShapeCoeffs,    // peak shape data
            double dmSearchRange, double dmFwhmRange    // needed for design matrix calculation
            )
        {
            Workspace W = new Workspace();

            W.SpectralData.MassAxis = massAxis;
            W.SpectralData.SignalAxisCrop = signalAxis;

            // TODO: fill the cluster dictionary


            if (resolutionPPBreaks == null)
            {
                W.Calibration.ResolutionInterp = new PolyInterpolation(resolutionFitCoeffs[0]);
            }
            else
            {
                W.Calibration.ResolutionInterp = new PPInterpolation(resolutionPPBreaks, resolutionFitCoeffs);
            }

            W.Calibration.Shape.Breaks = peakShapeBreaks;
            W.Calibration.Shape.Coeffs = peakShapeCoeffs;

            W.BuildDesignMatrix(dmSearchRange, dmFwhmRange);  // TODO: make the parameters adjustable

            W.FitAbundances();

            // TODO: extract the abundance (and errors maybe, anything else?) and put it to some dictionary, return that
        }

        public static void bla(int x)
        {
            Console.WriteLine(x);
        }

        public static void TestFit( /*Dictionary<string, Dictionary<string, object>> clusters*/ Dictionary<string, int> bla )
        {
            foreach (var item in bla.Values)
            {
                Console.WriteLine(item);
            }

            foreach (var key in bla.Keys)
            {
                Console.Write("Key: " + key);
                Console.WriteLine("\t\tVal: " + bla[key]);
            }

            //bla["jozko"] = 5;

            //Dictionary<string, Dictionary<string, object>> abd;
            //return abd;
        }

        public static void TestFit2(OrderedDictionary od)
        {
            try
            {
                foreach (var key in od.Keys)
                {
                    IFData.Cluster c = od[key] as IFData.Cluster;

                    Console.WriteLine(key);
                    Console.WriteLine(c.CentreOfMass);
                    Console.WriteLine(c.PeakData.Mass[0] + " " + c.PeakData.Abundance[0]);
                    Console.WriteLine(c.PeakData.Mass[1] + " " + c.PeakData.Abundance[1]);
                    Console.WriteLine(c.PeakData.Mass[2] + " " + c.PeakData.Abundance[2]);
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.Message);
            }
        }
    }
}
