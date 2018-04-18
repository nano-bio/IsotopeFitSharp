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
            OrderedDictionary clusters,    // cluster data
            double[] resolutionPPBreaks, double[][] resolutionCoeffs,    // resolution fit data
            double[] peakShapeBreaks, double[][] peakShapeCoeffs,    // peak shape data
            double dmSearchRange, double dmFwhmRange    // needed for design matrix calculation
            )
        {
            Workspace W = new Workspace();

            W.SpectralData.MassAxis = massAxis;
            W.SpectralData.SignalAxis = signalAxis;
            W.Clusters = clusters;

            if (resolutionPPBreaks == null)
            {
                W.Calibration.ResolutionInterp = new PolyInterpolation(resolutionCoeffs[0]);
            }
            else
            {
                W.Calibration.ResolutionInterp = new PPInterpolation(resolutionPPBreaks, resolutionCoeffs);
            }

            W.Calibration.Shape.Breaks = peakShapeBreaks;
            W.Calibration.Shape.Coeffs = peakShapeCoeffs;

            W.BuildDesignMatrix(dmSearchRange, dmFwhmRange);

            W.FitAbundances();
        }
    }
}
