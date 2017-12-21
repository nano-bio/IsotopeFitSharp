using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics;
using MathNet.Numerics.LinearRegression;

namespace IsotopeFit
{
    /// <summary>
    /// Derived class that handles polynomial interpolations and evaluation of data. 
    /// </summary>
    public class PolyInterpolation : Interpolation
    {
        #region Constructors

        /// <summary>
        /// Creates new interpolation object and calculates polynomial interpolation parameters.
        /// </summary>
        /// <param name="x">Array of x values.</param>
        /// <param name="y">Array of y values.</param>
        /// <param name="order">Order of polynomial.</param>
        public PolyInterpolation(double[] x, double[] y, int order)
        {
            xValues = x;
            yValues = y;

            if (order > 13)   //TODO: put this back to 10
            {
                throw new InterpolationException("Maximal order of polynomial interpolation is 10.");
            }

            Order = order;

            Calculate(x, y, order);
        }

        /// <summary>
        /// Creates new interpolation object and stores already known interpolation coeficients.
        /// </summary>
        /// <remarks>
        /// Input polynomial coefficients have to be sorted by increasing power from left to right in array.
        /// </remarks>
        /// <param name="coefs">Array of polynomial coefficients.</param>
        public PolyInterpolation(double[] coefs)
        {
            Coefs = coefs;
        }

        #endregion

        #region Properties

        public double[] Coefs { get; private set; }
        public int Order { get; private set; }

        #endregion

        #region Methods

        /// <summary>
        /// Uses polynomial least-square fit to calculate interpolation coeficients.
        /// </summary>
        /// <param name="x">Array of x values.</param>
        /// <param name="y">Array of y values.</param>
        /// <param name="order">Order of polynomial.</param>
        private void Calculate(double[] x, double[] y, int order)
        {
            double[] coefs = Fit.Polynomial(x, y, order, DirectRegressionMethod.QR);
            Coefs = coefs;
        }

        /// <summary>
        /// Initializes a polynomial evaluation of a selected point x.
        /// </summary>
        /// <param name="x">Single x value.</param>
        /// <returns>Single evaluated y value.</returns>
        public override double Evaluate(double x)
        {
            return PolyEvaluate(x);
        }

        /// <summary>
        /// Initializes a polynomial evaluation of a selected array of points x.
        /// </summary>
        /// <param name="x">Array of x values.</param>
        /// <returns>Array of evaluated y values.</returns>
        public override double[] Evaluate(double[] x)
        {
            double[] MultiEval = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                MultiEval[i] = PolyEvaluate(x[i]);
            }
            return MultiEval;
        }

        /// <summary>
        /// Evaluate a polynomial at point x.
        /// </summary>
        /// <param name="x">Single x value to be evaluated.</param>
        /// <returns>Single evaluated y value.</returns>
        private double PolyEvaluate(double x)
        {
            return MathNet.Numerics.Evaluate.Polynomial(x, Coefs);
        }

        #endregion
    }
}
