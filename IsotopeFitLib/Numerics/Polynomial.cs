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
        internal static Vector<double> PolynomialFit(double[] x, double[] y)
        {
            double[] coefs = Fit.Polynomial(x, y, 3, DirectRegressionMethod.QR);
            return Vector<double>.Build.DenseOfArray(coefs);
        }

        //TODO: matus
        internal static double PolynomialEval()
        {
            throw new NotImplementedException();
        }
    }
}
