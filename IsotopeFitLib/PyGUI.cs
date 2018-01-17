using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IsotopeFit
{
    public static class PyGUI
    {
        public static void Fit(
            double[] massAxis, double[] signalAxis,    // spectral data
            Dictionary<string, IFData.Cluster> clusters, List<string> clusterIDList,    // cluster data
            int resolutionFitType, double[] resolutionBreaks, double[] resolutionCoeffs,    // resolution fit data
            double[] peakShapeBreaks, double[] peakShapeCoeffs,    // peak shape data
            double dmSearchRange, double dmFwhmRange    // needed for design matrix calculation
            )
        {
            Workspace W = new Workspace();

            W.SpectralData.MassAxis = massAxis;
            W.SpectralData.SignalAxisCrop = signalAxis;
            
            // TODO: fill the cluster dictionary

            W.BuildDesignMatrix(dmSearchRange, dmFwhmRange);  // TODO: make the parameters adjustable

            W.FitAbundances();

            // TODO: extract the abundance (and errors maybe, anything else?) and put it to some dictionary, return that
        }
    }
}
