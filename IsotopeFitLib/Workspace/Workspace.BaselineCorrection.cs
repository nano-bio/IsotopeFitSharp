using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IsotopeFit
{
    /// <summary>
    /// This class is the main interface for outside usage.
    /// </summary>
    public partial class Workspace
    {
        /// <summary>
        /// Generates the x and y coordinates of the baseline break points.
        /// </summary>
        /// <param name="numOfSections">Optional, number of baseline break points to be generated.</param>
        /// <param name="cutoffLevel">Optional, baseline cutoff level specified as a double from 0 to 100. (Percentage of lowest y-values to be cut off.)</param>
        public void GenerateBaseline(int numOfSections = -1, double cutoffLevel = -1)
        {
            // update the baseline correction parameters, if they were supplied in the method call
            if (numOfSections != -1) BaselineCorrData.NumOfSections = numOfSections;
            if (cutoffLevel != -1) BaselineCorrData.CutoffLevel = cutoffLevel;

            // find the indices that indicate the interval marked by the supplied start and end mass
            int startIndex = Array.BinarySearch(SpectralData.RawMassAxis, BaselineCorrData.StartMass);
            if (startIndex < 0)
            {
                startIndex = ~startIndex;
            }

            int endIndex = Array.BinarySearch(SpectralData.RawMassAxis, BaselineCorrData.EndMass);
            if (endIndex < 0)
            {
                endIndex = ~endIndex - 1;   //because we want the last mass that is still smaller than the boundary
            }

            // cut the mass axis and signal
            double[] m = new double[endIndex - startIndex + 1];
            double[] y = new double[endIndex - startIndex + 1];

            Array.Copy(SpectralData.RawMassAxis, startIndex, m, 0, m.Length);
            Array.Copy(SpectralData.RawSignalAxis, startIndex, y, 0, y.Length);

            int step = (endIndex - startIndex) / BaselineCorrData.NumOfSections; //TODO: check if the index difference is equal to the array length
            int numOfValues = (int)Math.Floor(step * BaselineCorrData.CutoffLevel / 100);

            double[] corrXAxis = new double[BaselineCorrData.NumOfSections];
            double[] corrYAxis = new double[BaselineCorrData.NumOfSections];

            for (int i = 0; i < BaselineCorrData.NumOfSections; i++) //TODO: parallel for?
            {
                double[] s = y.Skip(i * step).Take(step).ToArray();
                Array.Sort(s);

                s = s.Take(numOfValues).ToArray();

                corrXAxis[i] = m[(2 * i + 1) * step / 2 + 1];  // the number will always be integer, even if there would be some decimal places in normal calculation
                corrYAxis[i] = s.Average();
            }

            BaselineCorrData.XAxis = corrXAxis;
            BaselineCorrData.YAxis = corrYAxis;
        }

        /// <summary>
        /// Calculates baseline corrected signal from raw signal data and baseline correction points. Stores the result in the <see cref="Workspace.SpectralData"/>.SignalAxis property.
        /// </summary>
        /// <remarks>
        /// <para>This method uses the PCHIP interpolation to obtain the baseline values, which are stored in the <see cref="Workspace.SpectralData"/>.Baseline property.</para>
        /// <para>Optional arguments <paramref name="xAxis"/> and <paramref name="yAxis"/> are meant to ease the process of loading data to the <see cref="Workspace"/>.
        /// If not supplied, the funcion will use previously stored values.</para>
        /// </remarks>
        /// <param name="xAxis">Optional x-axis for the baseline correction.</param>
        /// <param name="yAxis">Optional y-axis for the baseline correction.</param>
        /// <exception cref="WorkspaceException">Thrown when some of the required data have not been loaded in the <see cref="Workspace"/>.</exception>
        public void CorrectBaseline(double[] xAxis = null, double[] yAxis = null)
        {
            // raw spectral data are required
            if (SpectralData.RawSignalAxis == null || SpectralData.RawMassAxis == null) throw new WorkspaceException("Raw spectral data not specified.");
            if (SpectralData.RawSignalAxis.Length != SpectralData.RawMassAxis.Length) throw new WorkspaceException("Supplied spectral data axis have different lengths.");

            // if there are optional background correction points specified, save them to the workspace and use those
            if (xAxis != null) BaselineCorrData.XAxis = xAxis;
            if (yAxis != null) BaselineCorrData.YAxis = yAxis;

            // safety check, they can still be null in come use cases
            if (BaselineCorrData.XAxis == null || BaselineCorrData.YAxis == null) throw new WorkspaceException("Baseline correction points not specified.");
            if (BaselineCorrData.XAxis.Length != BaselineCorrData.YAxis.Length) throw new WorkspaceException("Supplied baseline correction point arrays have different lengths.");

            int massAxisLength = SpectralData.RawLength;

            //TODO: Evaluating the bg correction for the whole range might be useless. Specifiyng a mass range would make sense.
            //TODO: crop the mass axis to the last specified point of the baseline? might break usefulness.
            //PPInterpolation baselineFit = new PPInterpolation(BaselineCorrData.XAxis, BaselineCorrData.YAxis, PPInterpolation.PPType.PCHIP);    // in the matlab code it is also hard-coded pchip
            BaselineCorrData.BaselineInterpolation = new PPInterpolation(BaselineCorrData.XAxis, BaselineCorrData.YAxis, PPInterpolation.PPType.PCHIP);    // in the matlab code it is also hard-coded pchip

            //SpectralData.Baseline = new double[massAxisLength];
            SpectralData.SignalAxis = new double[massAxisLength];

            for (int i = 0; i < massAxisLength; i++)
            {
                //SpectralData.Baseline[i] = BaselineCorrData.BaselineInterpolation.Evaluate(SpectralData.RawMassAxis[i]);
                //SpectralData.SignalAxis[i] = SpectralData.RawSignalAxis[i] - SpectralData.Baseline[i];
                SpectralData.SignalAxis[i] = SpectralData.RawSignalAxis[i] - BaselineCorrData.BaselineInterpolation.Evaluate(SpectralData.RawMassAxis[i]);
            }
        }
    }
}
