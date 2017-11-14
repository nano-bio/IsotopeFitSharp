using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;

namespace IsotopeFit.Numerics
{
    public class Interpolation
    {
        #region Constructors

        public Interpolation(double[] x, double[] y, InterpolationType t)
        {
            Type = t;
            XValues = x;
            YValues = y;

            Calculated = Calculate(x, y, t);
        }

        public Interpolation(double[] coefs)
        {
            Type = InterpolationType.Polynomial;
            Calculated = true;
            Coefs = Matrix<double>.Build.DenseOfRowVectors(Vector<double>.Build.DenseOfArray(coefs));
        }

        public Interpolation(double[] breaks, double[][] coefs)
        {
            Type = InterpolationType.PiecewisePolynomial;
            Calculated = true;

            Breaks = Vector<double>.Build.DenseOfArray(breaks);
            Coefs = Matrix<double>.Build.DenseOfRowArrays(coefs);
        }

        #endregion

        public enum InterpolationType
        {
            Spline,
            PCHIP,
            Polynomial,
            PiecewisePolynomial
        }

        #region Properties

        public InterpolationType Type { get; private set; }
        public double[] XValues { get; private set; }
        public double[] YValues { get; private set; }
        public bool Calculated { get; private set; }
        public Matrix<double> Coefs { get; private set; }
        public Vector<double> Breaks { get; private set; }

        #endregion

        #region Methods

        private bool Calculate(double[] x, double[] y, InterpolationType t)
        {
            switch (t)
            {
                case InterpolationType.Spline:
                    break;
                case InterpolationType.PCHIP:
                    break;
                case InterpolationType.Polynomial:
                    break;
                default:
                    throw new InterpolationException("Unknown interpolation type.");
            }

            return true;
        }

        internal double Evaluate(double x)
        {
            if (!Calculated) throw new InterpolationException("Can not evaluate, interpolation object has not yet been supplied with input data. Use one of the constructors.");

            double retval;

            switch (Type)
            {
                case InterpolationType.Spline:
                case InterpolationType.PCHIP:
                case InterpolationType.PiecewisePolynomial:
                    //TODO: evaluate piecewise polynomial
                    break;
                case InterpolationType.Polynomial:
                    //TODO: evaluate polynomial
                    break;
                default:
                    throw new InterpolationException("Unknown interpolation type.");
            }

            return 0; //TODO
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
