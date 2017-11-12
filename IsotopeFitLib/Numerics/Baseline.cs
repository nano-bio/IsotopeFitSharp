﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Interpolation;

namespace IsotopeFit
{
    public static partial class Numerics
    {
        //TODO: methods for backround subtraction

        /// <summary>
        /// Corrects the raw data for baseline, according to supplied baseline correction information.
        /// </summary>
        /// <param name="bc">Data object containing the baseline correction information.</param>
        /// <param name="rd">Data object containing the raw experimental data.</param>
        /// <returns>MathNet vector containing the signal with the baseline subtracted from it.</returns>
        internal static Vector<double> CorrectBaseline(IFData.BaselineCorr bc, IFData.RawData rd)
        {
            int massAxisLength = rd.Length;

            //TODO: Evaluating the bg correction for the whole range might be useless. Specifiyng a mass range would make sense.
            Vector<double> baseline = Vector<double>.Build.DenseOfArray(SPPCHIP(bc.XAxis.ToArray(), bc.YAxis.ToArray(), rd.MassAxis.ToArray()));
            Vector<double> correctedSignal = Vector<double>.Build.Dense(massAxisLength, 0);

            for (int i = 0; i < massAxisLength; i++)
            {
                correctedSignal[i] = rd.SignalAxis[i] - baseline[i];
            }

            return correctedSignal;
        }

        /// <summary>
        /// Calculates shape preserving piecewise cubic hermite interpolation.
        /// </summary>
        /// <remarks>
        /// Before the interpolation itself, we need to calculate
        /// the derivative values to maintain interval monotonicity and function shape.
        /// Method used is from:
        /// F.N.Fritch, R.E.Carlson, 1980, Monotone Piecewise Cubic Interpolation, SIAM J. Numer.Anal. 17, pp. 238-246
        /// </remarks>
        /// <param name="x">Array of x values.</param>
        /// <param name="y">Array of y values.</param>
        /// <param name="xToEval">Array of x values, for which the interpolated curve is to be evaluated at.</param>
        /// <returns>Array of evaluated y values.</returns>
        private static double[] SPPCHIP(double[] x, double[] y, double[] xToEval)
        {
            /*
             * Before the interpolation itself, we need to calculate
             * the derivative values to maintain interval monotonicity and function shape.
             * Method used is from:
             * F. N. Fritch, R. E. Carlson, 1980, Monotone Piecewise Cubic Interpolation, SIAM J. Numer. Anal. 17, pp. 238-246
             * It is also used by GNU Octave.
             */

            // number of input and output values
            int n = x.Length;
            int nInterp = xToEval.Length;

            // variables related to horizontal data
            double h0, h1, hSum;

            // variables related to function deltas
            double delta0, delta1, deltaSave;
            double deltaMax, deltaMin;
            double deltaAdj0, deltaAdj1;

            // weights for calculating the derivative
            double w0, w1;

            double[] deriv = new double[n];
            
            // array for output values
            double[] interpVal = new double[nInterp];


            h0 = x[1] - x[0];
            delta0 = (y[1] - y[0]) / h0;
            deltaSave = delta0;

            h1 = x[2] - x[1];
            delta1 = (y[2] - y[1]) / h1;

            // 3-point formula
            hSum = h0 + h1;

            w0 = (h0 + hSum) / hSum;
            w1 = -h0 / hSum;

            deriv[0] = w0 * delta0 + w1 * delta1;

            if (Math.Sign(deriv[0]) * (Math.Sign(delta0)) <= 0)
            {
                deriv[0] = 0;
            }
            else if (Math.Sign(delta0) * Math.Sign(delta1) < 0)
            {
                deltaMax = 3.0 * delta0;
                if (Math.Abs(deriv[0]) > Math.Abs(deltaMax)) deriv[0] = deltaMax;
            }

            // internal points loop
            for (int i = 1; i < n - 1; i++)
            {
                if (i != 1)
                {
                    h0 = h1;
                    h1 = x[i + 1] - x[i];
                    hSum = h0 + h1;
                    delta0 = delta1;
                    delta1 = (y[i + 1] - y[i]) / h1;
                }

                deriv[i] = 0;
                if (Math.Sign(delta0) * Math.Sign(delta1) < 0)
                {
                    deltaSave = delta1;
                    continue;
                }

                if (Math.Sign(delta0) * Math.Sign(delta1) == 0)
                {
                    if (delta1 == 0) continue;
                    deltaSave = delta1;
                    continue;
                }

                w0 = (hSum + h0) / (3 * hSum);
                w1 = (hSum + h1) / (3 * hSum);

                deltaMax = Math.Max(Math.Abs(delta0), Math.Abs(delta1));
                deltaMin = Math.Min(Math.Abs(delta0), Math.Abs(delta1));

                deltaAdj0 = delta0 / deltaMax;
                deltaAdj1 = delta1 / deltaMax;

                deriv[i] = deltaMin / (w0 * deltaAdj0 + w1 * deltaAdj1);
            }

            // Special case - last point
            w0 = -h1 / hSum;
            w1 = (h1 + hSum) / hSum;
            deriv[n - 1] = w0 * delta0 + w1 * delta1;

            if (Math.Sign(deriv[n - 1]) * Math.Sign(delta1) <= 0)
            {
                deriv[n - 1] = 0;
            }
            else if (Math.Sign(delta0) * Math.Sign(delta1) < 0)
            {
                deltaMax = 3.0 * delta1;
                if (Math.Abs(deriv[n - 1]) > Math.Abs(deltaMax)) deriv[n - 1] = deltaMax;
            }
            
            // derivatives are calculated, let's interpolate
            CubicSpline cs = CubicSpline.InterpolateHermite(x, y, deriv);

            for (int i = 0; i < xToEval.Length; i++)
            {
                interpVal[i] = cs.Interpolate(xToEval[i]);
            }
          
            return interpVal;
        }
    }
}