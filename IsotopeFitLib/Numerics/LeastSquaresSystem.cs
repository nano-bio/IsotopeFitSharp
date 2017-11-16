﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;

namespace IsotopeFit.Numerics
{
    /// <summary>
    /// Class for handling least squares systems.
    /// </summary>
    public class LeastSquaresSystem
    {
        /// <summary>
        /// Creates new least squares system and populates it with the supplied data.
        /// </summary>
        /// <param name="desMat">Design matrix of the least squares system.</param>
        /// <param name="obsVec">Vector of observations of the least squares system.</param>
        public LeastSquaresSystem(Matrix<double> desMat, Vector<double> obsVec)
        {
            DesignMatrix = desMat;
            ObservationVector = obsVec;
        }

        internal Matrix<double> DesignMatrix { get; private set; }
        internal Vector<double> ObservationVector { get; private set; }
        public Vector<double> Solution { get; private set; }

        /// <summary>
        /// Calls the least square solver method and stores the result in the Solution property.
        /// </summary>
        /// <remarks>
        /// At the moment, only non-negative least squares solver is implemented.
        /// </remarks>
        public void Solve()
        {
            // at the moment we only need the NNLS method, so no need to add switches for more types
            Solution = NNLS(DesignMatrix, ObservationVector);
        }

        /// <summary>
        /// Calculate a non-negative least squares solution of C * x = d
        /// </summary>
        /// <remarks>
        /// The algorithm was originally published in:
        /// Lawson, Hanson, Solving Least Squares Problems, 1987, ISBN 978-0-89871-356-5, p. 160, Chapter 23.3
        /// Another description of the same algorithm, but easier to understand has been published in:
        /// Bro, De Jong, 1997, A fast non-negativity-constrained least squares algorithm, J. Chemom. 11, pp. 393-401
        /// </remarks>
        /// <param name="C">Matrix describing the model.</param>
        /// <param name="d">Vector with observation values.</param>
        /// <returns>MathNet vector with the solution.</returns>
        private static Vector<double> NNLS(Matrix<double> C, Vector<double> d)
        {
            //TODO: this needs to be set at the start of all calculations, right after Isotopefitter is called
            //MathNet.Numerics.Control.UseNativeMKL();
            //MathNet.Numerics.Control.UseMultiThreading();

            //Declarations and initializations, where necessary.
            int m = C.RowCount;
            int n = C.ColumnCount;

            // related to matrix C
            Matrix<double> CT = C.Transpose();
            double[][] Carr = C.ToColumnArrays();
            Matrix<double> CP;
            List<double[]> CPList = new List<double[]>();

            // solution vectors
            Vector<double> x = Vector<double>.Build.Dense(n, 0);
            Vector<double> z;

            // gradient vectors
            Vector<double> w;
            Vector<double> wz = Vector<double>.Build.Dense(n, 0);

            // active and passive sets
            bool[] P = new bool[n];
            bool[] A = new bool[n];

            for (int i = 0; i < n; i++)
            {
                P[i] = false;
                A[i] = true;
            }

            // helper variables
            Vector<double> resid;
            int outerIter = 0;
            int innerIter = 0;
            int iterMax = 3 * n;
            int zColIndex = 0;
            double tolx = 10 * MathNet.Numerics.Precision.DoublePrecision * n * C.L1Norm(); // that is how octave does it

            // calculation starts here

            System.Diagnostics.Debug.WriteLine(d.Count + " " + C.RowCount + " " + C.ColumnCount + " " + x.Count);

            resid = d - C * x;
            w = CT * resid;

            /* 
             * Main loop.
             * First condition checks whether there are still design parameters in the active set.
             * The second condition first filters out those elements from vector w, that correspond to active design parameters.
             * Then, those filtered values are compared to convergence tolerance.
             * In human language: In the end, the active set should be either empty, or the corresponding w elements must be close to zero.
             */
            while (A.Any(b => b == true) && w.Where((val, idx) => A[idx] == true).Any(val => val > tolx))
            {
                outerIter += 1;


                //TODO: check if wz initialization can be implemented using LINQ, something like wz.Where((val, idx) => P[idx] == false).Select(val => val = double.MinValue);
                for (int i = 0; i < n; i++)
                {
                    System.Diagnostics.Debug.Assert(P[i] != A[i]);

                    if (P[i])
                    {
                        wz[i] = double.MinValue;
                    }
                    else    // that means Z[i] is true
                    {
                        wz[i] = w[i];
                    }
                }

                // we move the design parameter with the highest gradient (w) to the passive set
                int t = wz.MaximumIndex();
                P[t] = true;
                A[t] = false;

                // find intermediary solution (z)
                CPList.Clear();

                for (int i = 0; i < n; i++)
                {
                    if (P[i])
                    {
                        CPList.Add(Carr[i]);
                    }
                }

                CP = Matrix<double>.Build.DenseOfColumnArrays(CPList.ToArray());
                z = CP.QR(QRMethod.Thin).Solve(d);

                //inner loop - check if any regression coefficient has turned negative
                while (z.Where((val, idx) => P[idx] == true).Any(val => val <= 0))
                {
                    innerIter++;
                    if (innerIter > iterMax)
                    {
                        Console.WriteLine("max pocet iteracii");
                    }

                    List<double> Q = new List<double>();

                    zColIndex = 0;

                    // get all negative regression coefficients
                    for (int i = 0; i < n; i++)
                    {
                        if (P[i])
                        {
                            if (z[zColIndex] <= 0)
                            {
                                //TODO: maybe it does not have to be a list, because we iterate through the inner loop
                                Q.Add(x[i] / (x[i] - z[zColIndex]));
                            }

                            zColIndex++;
                        }
                    }

                    // find the optimal alpha to make correction with
                    double alpha = Q.Min();
                    zColIndex = 0;

                    // find the solution that does not violate the constraints using that alpha
                    for (int i = 0; i < n; i++)
                    {
                        if (P[i])
                        {
                            x[i] += alpha * (z[zColIndex] - x[i]);
                            zColIndex++;
                        }
                        else
                        {
                            x[i] -= alpha * x[i];
                        }

                        // update P and Z accordingly
                        if (Math.Abs(x[i]) < tolx && P[i])
                        {
                            A[i] = true;
                        }

                        P[i] = !A[i];
                    }

                    // recalculate for z
                    CPList.Clear();

                    for (int i = 0; i < n; i++)
                    {
                        if (P[i])
                        {
                            CPList.Add(Carr[i]);
                        }
                    }

                    CP = Matrix<double>.Build.DenseOfColumnArrays(CPList.ToArray());
                    z = CP.QR(QRMethod.Thin).Solve(d);
                }

                // calculate gradient
                zColIndex = 0;

                for (int i = 0; i < n; i++)
                {
                    if (P[i])
                    {
                        x[i] = z[zColIndex];
                        zColIndex++;
                    }
                }

                resid = d - C * x;
                w = CT * resid;
            }

            return x;
        }
    }
}
