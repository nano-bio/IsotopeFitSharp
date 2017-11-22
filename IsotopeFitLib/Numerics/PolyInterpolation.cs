using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics;
using MathNet.Numerics.LinearRegression;

namespace IsotopeFit.Numerics
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
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="order"></param>
        public PolyInterpolation(double[] x, double[] y, int order)
        {
            xValues = x;
            yValues = y;

            if (order>10)
            {
                throw new InterpolationException("Maximal order of polynomial interpolation is 10.");
            }

            Order = order;

            Calculate(x, y, order);
        }
        
        /// <summary>
        /// Creates new interpolation object and stores already known interpolation coeficients.
        /// </summary>
        /// <param name="coefs"></param>
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
        /// Uses polynomial least-square fit to calculate interpolation coeficients 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="order"></param>
        private void Calculate(double[] x, double[] y, int order)
        {
            double[] coefs = Fit.Polynomial(x, y, order, DirectRegressionMethod.QR);
            Coefs = coefs;
        }
        
        public override double Evaluate(double x)
        {
            return PolyEvaluate(x);
        }

        public override double[] Evaluate(double[] x)
        {
            double[] MultiEval = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                MultiEval[i] = PolyEvaluate(x[i]);
            }
            return MultiEval;
        }

        private double PolyEvaluate(double x)
        {
            return MathNet.Numerics.Evaluate.Polynomial(x, Coefs);
        }

        #endregion
    }
}
