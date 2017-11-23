using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Interpolation;
using MathNet.Numerics.Properties;
using System.IO;

namespace IsotopeFit.Numerics
{
    /// <summary>
    /// Derived class that handles piecewise polynomial interpolations and evaluation of data. 
    /// </summary>
    public class PPInterpolation : Interpolation
    {
        #region Constructors

        /// <summary>
        /// Creates new interpolation object and calculates piecewise polynomial interpolation parameters.
        /// </summary>
        /// <param name="x">Array of x values.</param>
        /// <param name="y">Array of y values.</param>
        /// <param name="t">Piecewise polynomial interpolation type.</param>
        public PPInterpolation(double[] x, double[] y, PPType t)
        {
            xValues = x;
            yValues = y;
            ppType = t;

            Calculate(x, y, t);
        }

        /// <summary>
        /// Creates new interpolation object and stores already known interpolation parameters.
        /// </summary>
        /// <remarks>
        /// Input polynomial coefficients have to be sorted by increasing power from left to right in 2D array.
        /// </remarks>
        /// <param name="breaks">Array of break points.</param>
        /// <param name="coefs">2D array of polynomial coefficients.</param>
        public PPInterpolation(double[] breaks, double[][] coefs)
        {
            Breaks = breaks;
            Coefs = coefs;
        }

        #endregion

        #region Properties

        public PPType ppType { get; private set; }
        public double[] Breaks { get; private set; }
        public double[][] Coefs { get; private set; }

        #endregion

        public enum PPType
        {
            Spline,
            PCHIP
        }

        #region Methods

        /// <summary>
        /// Initializes a calculation by selected piecewise polynomial interpolation method. 
        /// </summary>
        /// <param name="x">Array of x values.</param>
        /// <param name="y">Array of y values.</param>
        /// <param name="t">Piecewise polynomial interpolation type.</param>
        private void Calculate(double[] x, double[] y, PPType t)
        {
            switch (t)
            {
                case PPType.Spline:
                    Spline(x, y);
                    break;
                case PPType.PCHIP:
                    PCHIP(x, y);
                    break;
                default:
                    throw new InterpolationException("Unknow interpolation type.");
            }
        }
        
        //TODO: implement spline calculation - continuous 2nd derivation
        internal static Matrix<double> Spline(double[] x, double[] y)
        {
            throw new NotImplementedException();
        }
        
        /// <summary>
        /// Calculates shape preserving piecewise cubic hermite polynomial interpolation.
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
        private void PCHIP(double[] x, double[] y)
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
            //int nInterp = xToEval.Length;

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
            //double[] interpVal = new double[nInterp];


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
            
            // Create a hermite cubic spline interpolation from a set of (x,y) value pairs and their slope (first derivative), sorted ascendingly by x.
            if (x.Length != y.Length || x.Length != deriv.Length)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            if (x.Length < 2)
            {
                throw new ArgumentException(string.Format(Resources.ArrayTooSmall, 2), "x");
            }

            // Constant, liner, quadratic and cubic coefficient, respectively.
            double[] c0 = new double[x.Length - 1];
            double[] c1 = new double[x.Length - 1];
            double[] c2 = new double[x.Length - 1];
            double[] c3 = new double[x.Length - 1];

            for (int i = 0; i < c0.Length; i++)
            {
                double w = x[i + 1] - x[i];
                double w2 = w * w;
                c0[i] = y[i];
                c1[i] = deriv[i];
                c2[i] = (3 * (y[i + 1] - y[i]) / w - 2 * deriv[i] - deriv[i + 1]) / w;
                c3[i] = (2 * (y[i] - y[i + 1]) / w + deriv[i] + deriv[i + 1]) / w2;
            }

            Breaks = x;
            Coefs = Matrix<double>.Build.DenseOfColumnArrays(c0, c1, c2, c3).ToRowArrays();
        }

        /// <summary>
        /// Initializes a piecewise cubic polynomial evaluation of a selected point x.
        /// </summary>
        /// <param name="x">Single x value.</param>
        /// <returns>Single evaluated y value.</returns>
        public override double Evaluate(double x)
        {
            return PPEval(Breaks, Coefs, x);
        }

        /// <summary>
        /// Initializes a piecewise cubic polynomial evaluation of a selected array of points x.
        /// </summary>
        /// <param name="x">Array of x values.</param>
        /// <returns>Array of evaluated y values.</returns>
        public override double[] Evaluate(double[] x)
        {
            double[] multiEval = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                multiEval[i] = PPEval(Breaks, Coefs, x[i]);
            }
            return multiEval;
        }

        /// <summary>
        /// Evaluates a piecewise cubic polynomial at point x within the defined range. 
        /// </summary>
        /// <remarks>
        /// Evaluation itself does not depend on selected calculation method.
        /// </remarks>
        /// <param name="breaks">Array of break points.</param>
        /// <param name="coefs">2D array of polynomial coefficients.</param>
        /// <param name="x">Single x value to be evaluated.</param>
        /// <returns>Single evaluated y value.</returns>
        internal static double PPEval(double[] breaks, double[][] coefs, double x)
        {
            // Assigns an index of equal break point or a bitwise compelent index of first greater one to the point x. 
            int breakIndex = Array.BinarySearch(breaks, x);

            // Selected point equals to one of the break points.
            if (breakIndex >= 0)
            {
                if (breakIndex == breaks.Length - 1)
                {
                    return MathNet.Numerics.Evaluate.Polynomial((x - breaks[breakIndex - 1]), coefs[breakIndex - 1]);
                }
                return coefs[breakIndex][0];
            }
            // Selected point is in between of two break points. 
            else
            {
                int indexOfGreater = ~breakIndex;

                if (indexOfGreater == 0) indexOfGreater++;
                if (indexOfGreater == breaks.Length) indexOfGreater--;
                
                if ((0 < indexOfGreater) && (indexOfGreater < breaks.Length))
                {
                    return MathNet.Numerics.Evaluate.Polynomial((x - breaks[indexOfGreater - 1]), coefs[indexOfGreater - 1]);
                }
                else
                {
                    throw new InterpolationException("Evaluation point is outside of definition range.");
                }
            }
            
            throw new NotImplementedException(); // TODO what is a purpose of the line?
        }

        #endregion
    }
}
