using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Interpolation;

namespace IsotopeFit.Numerics
{
    public static partial class Algorithm
    {
        //TODO: move to workspace

        /// <summary>
        /// Corrects the raw data for baseline, according to supplied baseline correction information.
        /// </summary>
        /// <param name="bc">Data object containing the baseline correction information.</param>
        /// <param name="rd">Data object containing the raw experimental data.</param>
        /// <returns>MathNet vector containing the signal with the baseline subtracted from it.</returns>
        internal static Vector<double> CorrectBaseline(IFData.BaselineCorr bc, IFData.Spectrum rd)
        {
            int massAxisLength = rd.Length;

            //TODO: Evaluating the bg correction for the whole range might be useless. Specifiyng a mass range would make sense.
            
            /*
             * Now we have to first create the interpolation object, which also calculates the interpolation parameters.
             * Only after that can we call the interpolation evaluation function.
             */
            
            //TODO: rewrite for the new interpolation scheme

            //Vector<double> baseline = Vector<double>.Build.DenseOfArray(PCHIP(bc.XAxis.ToArray(), bc.YAxis.ToArray(), rd.MassAxis.ToArray()));
            Vector<double> correctedSignal = Vector<double>.Build.Dense(massAxisLength, 0);

            for (int i = 0; i < massAxisLength; i++)
            {
                //correctedSignal[i] = rd.SignalAxis[i] - baseline[i];
            }

            return correctedSignal;
        }        
    }
}
