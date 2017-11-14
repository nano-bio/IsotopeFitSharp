using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearRegression;

namespace IsotopeFit.Numerics
{
    public static partial class Algorithm
    {
        internal static Vector<double> PolynomialFit(Vector<double> x, Vector<double> y)
        {
            double[] coefs = Fit.Polynomial(x.ToArray(), y.ToArray(), 3, DirectRegressionMethod.QR);
            return Vector<double>.Build.DenseOfArray(coefs);
        }

        //TODO: matus
        internal static double PolynomialEval()
        {
            throw new NotImplementedException();
        }
    }
}
