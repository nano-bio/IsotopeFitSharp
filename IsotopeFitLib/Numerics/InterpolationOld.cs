using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;

namespace IsotopeFit.Numerics
{
    /// <summary>
    /// Class that handles data interpolations and related data.
    /// </summary>
    public class InterpolationOld
    {
        #region Constructors

        /// <summary>
        /// Creates new interpolation object and calculates interpolation parameters.
        /// </summary>
        /// <param name="x">Array of sorted x-values that are to be interpolated.</param>
        /// <param name="y">Array of y-values that are to be interpolated.</param>
        /// <param name="t">Interpolation type.</param>
        public InterpolationOld(double[] x, double[] y, InterpType t)
        {
            Type = t;
            XValues = x;
            YValues = y;

            if (t != InterpType.Polynomial)
            {
                Breaks = Vector<double>.Build.DenseOfArray(x);
            }

            Calculated = Calculate(x, y, t);
        }

        /// <summary>
        /// Creates new polynomial interpolation object from polynomial coefficients.
        /// </summary>
        /// <param name="coefs"></param>
        public InterpolationOld(double[] coefs)
        {
            Type = InterpType.Polynomial;
            Calculated = true;
            Coefs = Matrix<double>.Build.DenseOfRowVectors(Vector<double>.Build.DenseOfArray(coefs));
        }

        public InterpolationOld(double[] breaks, double[][] coefs)
        {
            Type = InterpType.PiecewisePolynomial;
            Calculated = true;

            Breaks = Vector<double>.Build.DenseOfArray(breaks);
            Coefs = Matrix<double>.Build.DenseOfRowArrays(coefs);
        }

        #endregion

        public enum InterpType
        {
            Spline,
            PCHIP,
            Polynomial,
            PiecewisePolynomial
        }

        #region Properties

        public InterpType Type { get; private set; }
        public double[] XValues { get; private set; }
        public double[] YValues { get; private set; }
        public bool Calculated { get; private set; }
        public Matrix<double> Coefs { get; private set; }
        public Vector<double> Breaks { get; private set; }

        #endregion

        #region Methods

        private bool Calculate(double[] x, double[] y, InterpType t)
        {
            switch (t)
            {
                case InterpType.Spline:
                    Coefs = Algorithm.Spline(x, y);
                    break;
                case InterpType.PCHIP:
                    Coefs = Algorithm.PCHIP(x, y);
                    break;
                case InterpType.Polynomial:
//                    Coefs = Matrix<double>.Build.DenseOfRowVectors(Algorithm.PolynomialFit(x, y));
                    break;
                default:
                    throw new InterpolationException("Unknown interpolation type.");
            }

            return true;
        }

        internal double Evaluate(double x)
        {
            if (!Calculated) throw new InterpolationException("Can not evaluate, interpolation object has not yet been supplied with input data. Use one of the constructors.");

            switch (Type)
            {
                case InterpType.Spline:
                case InterpType.PCHIP:
                case InterpType.PiecewisePolynomial:
                    return Algorithm.PPEval(Breaks, Coefs, x);
 //               case InterpType.Polynomial:
 //                   return Algorithm.PolynomialEval();
                default:
                    throw new InterpolationException("Unknown interpolation type.");
            }
        }

        [Obsolete]
        internal double[] Evaluate(double[] x)
        {
            if (!Calculated) throw new InterpolationException("Can not evaluate, interpolation object has not yet been supplied with input data. Use one of the constructors.");

            double[] eval = new double[x.Length];

            switch (Type)
            {
                case InterpType.Spline:
                case InterpType.PCHIP:
                case InterpType.PiecewisePolynomial:
                    {
                        for (int i = 0; i < x.Length; i++)      //TODO: replace these for loops by implementing evaluation functions that take arrays as inputs - less function calls
                        {
                            eval[i] = Algorithm.PPEval(Breaks, Coefs, x[i]);
                        }
                        break;
                    }                    
                case InterpType.Polynomial:
                    {
                        for (int i = 0; i < x.Length; i++)
                        {
 //                           eval[i] = Algorithm.PolynomialEval();
                        }
                        break;
                    }
                default:
                    throw new InterpolationException("Unknown interpolation type.");
            }

            return eval;
        }

        #endregion

        [Serializable]
        public class InterpolationException : Exception
        {
            public InterpolationException() { }
            public InterpolationException(string message) : base(message) { }
            public InterpolationException(string message, Exception inner) : base(message, inner) { }
            protected InterpolationException(
              System.Runtime.Serialization.SerializationInfo info,
              System.Runtime.Serialization.StreamingContext context) : base(info, context) { }
        }
    }
}
