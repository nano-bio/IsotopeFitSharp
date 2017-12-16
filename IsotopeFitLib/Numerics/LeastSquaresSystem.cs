﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;

using CSparse.Double;

namespace IsotopeFit
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
        public LeastSquaresSystem(SparseMatrix desMat, Vector<double> obsVec)
        {
            DesignMatrix = desMat;
            ObservationVector = obsVec;
        }

        internal SparseMatrix DesignMatrix { get; private set; }
        internal Vector<double> ObservationVector { get; private set; }
        internal int[] ColumnOrdering { get; set; }
        public Vector<double> Solution { get; private set; }

        public void Solve2()
        {
            // partitioning 1: by non-overlaping diagonal elements
            // partitioning 2: by non-overlapping blocks(there does not have to be an element from case 1 between them)
            // partitioning 3: overlapping blocks cut duplicitly, some columns get fitted in both blocks (maybe three: left-center-right) - needs to be tested if correct

            double[] values = DesignMatrix.Values;
            int[] rowIndices = DesignMatrix.RowIndices;
            int[] columnPointers = DesignMatrix.ColumnPointers;

            SparseMatrix AT = (SparseMatrix)DesignMatrix.Transpose();

            double[] valuesT = AT.Values;
            int[] columnIndicesT = AT.RowIndices;
            int[] rowPointersT = AT.ColumnPointers;

            double[] colCounts = new double[columnPointers.Length - 1];
            double[] rowCounts = new double[rowPointersT.Length - 1];

            for (int i = 0; i < colCounts.Length; i++)
            {
                colCounts[i] = columnPointers[i + 1] - columnPointers[i];
                rowCounts[i] = rowPointersT[i + 1] - rowPointersT[i];
            }

            List<int> nonOverlap = new List<int>();
            List<int> cutCoordinates = new List<int>(); 

            for (int i = 0; i < colCounts.Length; i++)
            {
                if (rowCounts[i] == 1 && colCounts[i] == 1)     // partitioning 1
                {
                    nonOverlap.Add(i);
                    cutCoordinates.Add(i);
                }
                else if (i < colCounts.Length - 1 && rowCounts[i] == 1 && colCounts[i + 1] == 1)    // partitioning 2, it must not be the last row of the matrix
                {
                    cutCoordinates.Add(i);
                }
                else if (i == colCounts.Length - 1)    // if it is the last row (and it didnt fall into the first case), end the block
                {
                    cutCoordinates.Add(i);
                }
            }

            // now we check the resulting block sizes and if some of them are too big, we will apply the type 3 partitioning - but that sacrifices precision
            int[] blSizes = new int[cutCoordinates.Count];
            blSizes[0] = cutCoordinates[0] + 1;

            for (int i = 1; i < blSizes.Length; i++)
            {
                blSizes[i] = cutCoordinates[i] - cutCoordinates[i - 1]; //TODO: currently it is block sizes and blocks ends - maybe block beginnins would be nicer?
            }

            //int maxSize = blSizes.Max();

            // calculate the non-overlaping cluster abundances
            double[] abd = new double[DesignMatrix.ColumnCount];

            //foreach (int i in nonOverlap)
            //{
            //    abd[i] = ObservationVector[i] / DesignMatrix.At(i, i);
            //}

            for (int i = 0; i < cutCoordinates.Count; i++)
            {
                if (blSizes[i] == 1)
                {
                    abd[cutCoordinates[i]] = ObservationVector[cutCoordinates[i]] / DesignMatrix.At(cutCoordinates[i], cutCoordinates[i]);
                }
                else    // create the diagonal block and use NNLS
                {
                    int[] indexes = Enumerable.Range(cutCoordinates[i] - blSizes[i] + 1, blSizes[i]).ToArray();
                    int start = cutCoordinates[i] - blSizes[i] + 1;
                    int end = start + blSizes[i];

                    int valStart = columnPointers[start];
                    int valLength = columnPointers[end] - valStart;

                    double[] vals = new double[valLength];
                    int[] rIdx = new int[valLength];
                    int[] cPt = new int[end - start + 1];

                    Array.Copy(values, valStart, vals, 0, valLength);
                    Array.Copy(rowIndices, valStart, rIdx, 0, valLength);
                    Array.Copy(columnPointers, start, cPt, 0, end - start + 1);

                    int ggg = rIdx[0];
                    for (int j = 0; j < rIdx.Length; j++)
                    {
                        rIdx[j] -= ggg;
                    }

                    int hhh = cPt[0];
                    for (int j = 0; j < cPt.Length; j++)
                    {
                        cPt[j] -= hhh;
                    }

                    SparseMatrix C = new SparseMatrix(cPt.Length - 1, cPt.Length - 1);
                    {
                        C.Values = vals;
                        C.RowIndices = rIdx;
                        C.ColumnPointers = cPt;
                    };

                    double[] obs = new double[cPt.Length - 1];
                    Array.Copy(ObservationVector.ToArray(), start, obs, 0, end - start);

                    try
                    {
                        double[] sol = NNLS3(C, Vector<double>.Build.DenseOfArray(obs)).ToArray();
                        Array.Copy(sol, 0, abd, start, sol.Length);
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine(ex.Message);
                    }                    
                }
            }

            Solution = Vector<double>.Build.DenseOfArray(abd);

            Solution = Vector<double>.Build.DenseOfArray(Enumerable.Zip(ColumnOrdering, Solution, (idx, val) => new { idx, val }).OrderBy(v => v.idx).Select(v => v.val).ToArray());
        }

        /// <summary>
        /// Calls the least square solver method and stores the result in the Solution property.
        /// </summary>
        /// <remarks>
        /// At the moment, only non-negative least squares solver is implemented.
        /// </remarks>
        [Obsolete]
        public void Solve()
        {
            // at the moment we only need the NNLS method, so no need to add switches for more types
            Solution = NNLS3(DesignMatrix, ObservationVector);

            // reorder the solution according to the ColumnOrdering information from sparse QR factorization
            Solution = Vector<double>.Build.DenseOfArray(Enumerable.Zip(ColumnOrdering, Solution, (idx, val) => new { idx, val }).OrderBy(v => v.idx).Select(v => v.val).ToArray());
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
        private static Vector<double> NNLS3(SparseMatrix C, Vector<double> d)
        {
            //TODO: this function can be rewritten to work with sparse matrices and vectors - faster, less ram consumption, nicer

            //TODO: this can be set here as well, we dont really call that many mathnet functions
            //MathNet.Numerics.Control.UseNativeMKL();
            //MathNet.Numerics.Control.UseMultiThreading();

            //Declarations and initializations, where necessary.
            int m = C.RowCount;
            int n = C.ColumnCount;

            // related to matrix C
            SparseMatrix CT = (SparseMatrix)C.Transpose();
            double[][] Carr = new double[C.ColumnCount][];  //TODO: not effective


            for (int i = 0; i < C.ColumnCount; i++)
            {
                Carr[i] = C.Column(i);
            }

            SparseMatrix CP;
            List<double[]> CPList = new List<double[]>();
            //SparseMatrix CPs;
            double[][] CPArr = new double[C.ColumnCount][];

            //for (int i = 0; i < C.ColumnCount; i++)
            //{
            //    CPArr[i] = new double[C.ColumnCount];
            //}

            //CP = Matrix<double>.Build.DenseOfColumnArrays(CPArr);

            // solution vectors
            double[] x = new double[n];
            //double[] x = new double[n];
            Vector<double> temp = Vector<double>.Build.Dense(n, 0);
            double[] z = new double[n];

            // gradient vectors
            Vector<double> w = Vector<double>.Build.Dense(n, 0);
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

            //System.Diagnostics.Debug.WriteLine(d.Count + " " + C.RowCount + " " + C.ColumnCount + " " + x.Count);

            C.Multiply(x, temp.AsArray());  //TODO: may return null
            resid = d - temp;
            //w = CT * resid;
            CT.Multiply(resid.AsArray(), w.AsArray());  //TODO: may return null

            //Pardiso.Initialize(C);

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
                Array.Clear(z, 0, z.Length);

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
                //// count how many active columns we have
                //int pCols = 0;
                //for (int i = 0; i < n; i++)
                //{
                //    if (P[i]) pCols++;
                //}

                //// construct the colPointers array
                //int[] AcolPointers = new int[pCols + 1];

                //pCols = 1;
                //for (int i = 0; i < n; i++)
                //{
                //    if (P[i])
                //    {
                //        AcolPointers[pCols] = AcolPointers[pCols - 1] + C.ColumnPointers[i + 1] - C.ColumnPointers[i];
                //        pCols++;
                //    }
                //}

                //// now we can construct the values and row indices
                //double[] Avalues = new double[AcolPointers.Last()];
                //int[] ArowIndices = new int[AcolPointers.Last()];

                //pCols = 0;
                //for (int i = 0; i < n; i++)
                //{
                //    if (P[i])
                //    {
                //        int sectionLength = AcolPointers[pCols + 1] - AcolPointers[pCols];

                //        Array.Copy(C.Values, C.ColumnPointers[i], Avalues, AcolPointers[pCols], sectionLength);
                //        Array.Copy(C.RowIndices, C.ColumnPointers[i], ArowIndices, AcolPointers[pCols], sectionLength);

                //        pCols++;
                //    }
                //}

                //CP = new SparseMatrix(n, pCols)
                //{
                //    Values = Avalues,
                //    RowIndices = ArowIndices,
                //    ColumnPointers = AcolPointers
                //};

                CP = BuildCP(C, P);

                double[] zk = SPQR.Solve(CP, d.ToArray());

                zColIndex = 0;
                for (int i = 0; i < n; i++)
                {
                    if (P[i])
                    {
                        z[i] = zk[zColIndex];
                        zColIndex++;
                    }
                }

                //inner loop - check if any regression coefficient has turned negative
                while (z.Where((val, idx) => P[idx] == true).Any(val => val <= 0))
                {
                    innerIter++;
                    if (innerIter > iterMax)
                    {
                        throw new Exception("max pocet iteracii");
                    }

                    List<double> Q = new List<double>();

                    zColIndex = 0;

                    // get all negative regression coefficients
                    for (int i = 0; i < n; i++)
                    {
                        if (P[i] && z[i] <= 0)  // 
                        {
                            //if (z[zColIndex] <= 0) //
                            //{
                            //TODO: maybe it does not have to be a list, because we iterate through the inner loop
                            double derp = x[i] / (x[i] - z[i]);     //zColIndex
                            Q.Add(derp);
                            //}

                            //zColIndex++;
                        }
                    }

                    // find the optimal alpha to make correction with
                    double alpha = Q.Min();
                    //zColIndex = 0;

                    // find the solution that does not violate the constraints using that alpha
                    for (int i = 0; i < n; i++)
                    {
                        //if (P[i])
                        //{
                        x[i] += alpha * (z[i] - x[i]); //zColIndex
                                                       //zColIndex++;
                                                       //}
                                                       //else
                                                       //{
                                                       //    x[i] -= alpha * x[i];
                                                       //}

                        // update P and Z accordingly
                        if (Math.Abs(x[i]) < tolx && P[i])
                        {
                            A[i] = true;
                        }

                        P[i] = !A[i];
                    }

                    // recalculate for z
                    //CPList.Clear();
                    //Array.Copy(Carr, CPArr, Carr.Length);

                    //for (int i = 0; i < n; i++)
                    //{
                    //    if (P[i])
                    //    {
                    //        CPList.Add(Carr[i]);
                    //    }

                    //    //if (!P[i])
                    //    //{
                    //    //    for (int j = 0; j < CPArr[i].Length; j++)
                    //    //    {
                    //    //        CPArr[i][j] = 0;
                    //    //    }
                    //    //}
                    //}

                    //CP = Matrix<double>.Build.DenseOfColumnArrays(CPList.ToArray());
                    //CP = Matrix<double>.Build.DenseOfColumnArrays(CPArr);
                    //CP = Matrix<double>.Build.SparseOfColumnArrays(CPList.ToArray());
                    //zk = CP.QR(QRMethod.Thin).Solve(d).ToArray();
                    //zk = Pardiso.Solve(CP, d.ToArray(), P);
                    //z = Vector<double>.Build.Dense(ggg);

                    CP = BuildCP(C, P);

                    zk = SPQR.Solve(CP, d.ToArray());

                    Array.Clear(z, 0, z.Length);
                    zColIndex = 0;
                    for (int i = 0; i < n; i++)
                    {
                        if (P[i])  // && z[i] <= 0
                        {
                            z[i] = zk[zColIndex];   //
                            zColIndex++;
                        }
                    }

                    //for (int i = 0; i < z.Count; i++)
                    //{
                    //    if (!P[i]) z[i] = 0;
                    //}
                }

                // calculate gradient
                //zColIndex = 0;

                //for (int i = 0; i < n; i++)
                //{
                //if (P[i])
                //{
                //x[i] = z[i];    //zColIndex
                //zColIndex++;
                //}
                //}
                Array.Copy(z, x, z.Length);
                //x = Vector<double>.Build.DenseOfVector(z);

                //resid = d - C * x;
                //w = CT * resid;

                try
                {
                    C.Multiply(x, temp.AsArray());  //TODO: may return null
                    resid = d - temp;
                    //w = CT * resid;
                    CT.Multiply(resid.AsArray(), w.AsArray());  //TODO: may return null
                }
                catch (Exception)
                {

                    throw;
                }
            }

            return Vector<double>.Build.DenseOfArray(x);
        }

        private static SparseMatrix BuildCP(SparseMatrix C, bool[] P)
        {
            // count how many active columns we have
            int pCols = 0;
            for (int i = 0; i < C.RowCount; i++)
            {
                if (P[i]) pCols++;
            }

            // construct the colPointers array
            int[] AcolPointers = new int[pCols + 1];

            pCols = 1;
            for (int i = 0; i < C.RowCount; i++)
            {
                if (P[i])
                {
                    AcolPointers[pCols] = AcolPointers[pCols - 1] + C.ColumnPointers[i + 1] - C.ColumnPointers[i];
                    pCols++;
                }
            }

            // now we can construct the values and row indices
            double[] Avalues = new double[AcolPointers.Last()];
            int[] ArowIndices = new int[AcolPointers.Last()];

            pCols = 0;
            for (int i = 0; i < C.RowCount; i++)
            {
                if (P[i])
                {
                    int sectionLength = AcolPointers[pCols + 1] - AcolPointers[pCols];

                    Array.Copy(C.Values, C.ColumnPointers[i], Avalues, AcolPointers[pCols], sectionLength);
                    Array.Copy(C.RowIndices, C.ColumnPointers[i], ArowIndices, AcolPointers[pCols], sectionLength);

                    pCols++;
                }
            }

            SparseMatrix CP = new SparseMatrix(C.RowCount, pCols)
            {
                Values = Avalues,
                RowIndices = ArowIndices,
                ColumnPointers = AcolPointers
            };

            return CP;
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
        [Obsolete]
        private static Vector<double> NNLS2(SparseMatrix C, Vector<double> d)
        {
            //TODO: this function can be rewritten to work with sparse matrices and vectors - faster, less ram consumption, nicer

            //TODO: this can be set here as well, we dont really call that many mathnet functions
            //MathNet.Numerics.Control.UseNativeMKL();
            //MathNet.Numerics.Control.UseMultiThreading();

            //Declarations and initializations, where necessary.
            int m = C.RowCount;
            int n = C.ColumnCount;

            // related to matrix C
            SparseMatrix CT = (SparseMatrix)C.Transpose();
            double[][] Carr = new double[C.ColumnCount][];  //TODO: not effective


            for (int i = 0; i < C.ColumnCount; i++)
            {
                Carr[i] = C.Column(i);
            }

            Matrix<double> CP;
            List<double[]> CPList = new List<double[]>();
            //SparseMatrix CPs;
            double[][] CPArr = new double[C.ColumnCount][];

            //for (int i = 0; i < C.ColumnCount; i++)
            //{
            //    CPArr[i] = new double[C.ColumnCount];
            //}

            //CP = Matrix<double>.Build.DenseOfColumnArrays(CPArr);

            // solution vectors
            double[] x = new double[n];
            //double[] x = new double[n];
            Vector<double> temp = Vector<double>.Build.Dense(n, 0);
            double[] z = new double[n];

            // gradient vectors
            Vector<double> w = Vector<double>.Build.Dense(n, 0);
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

            //System.Diagnostics.Debug.WriteLine(d.Count + " " + C.RowCount + " " + C.ColumnCount + " " + x.Count);

            C.Multiply(x, temp.AsArray());  //TODO: may return null
            resid = d - temp;
            //w = CT * resid;
            CT.Multiply(resid.AsArray(), w.AsArray());  //TODO: may return null

            //Pardiso.Initialize(C);

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
                Array.Clear(z, 0, z.Length);

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

                // pre dlzku colPointers musime najprv zratat kolko je pasivnych stlpcov
                //int pCols = 0;
                //for (int i = 0; i < n; i++)
                //{
                //    if (P[i]) pCols++;
                //}

                // kolko bude nenulovych elementov? no nie viac ako povodny pocet nenulovych

                //double[] values = new double[C.NonZerosCount];
                //int[] rowIndices = new int[C.NonZerosCount];
                //int[] colPointers = new int[pCols];

                //colPointers[0] = 0;

                //Array.Copy(Carr, CPArr, Carr.Length);

                for (int i = 0; i < n; i++)
                {
                    //int numOfValsToRead = C.ColumnPointers[i + 1] - C.ColumnPointers[i];
                    //int columnBeginIndex = C.ColumnPointers[i];

                    if (P[i])
                    {
                        CPList.Add(Carr[i]);

                        //Array.Copy(C.Values, columnBeginIndex, values, columnBeginIndex, numOfValsToRead);
                        //Array.Copy(C.RowIndices, columnBeginIndex, rowIndices, columnBeginIndex, numOfValsToRead);
                    }

                    //if (!P[i])
                    //{
                    //    for (int j = 0; j < CPArr[i].Length; j++)
                    //    {
                    //        CPArr[i][j] = 0;
                    //    }
                    //}

                    //else    // set the rest to zero, because the columns might have been reactivated
                    //{
                    //    // TODO: zatial netreba, lebo teraz sa tie arraye robia vzdy nanovo, potom ich presuniem mimo slucku
                    //    //for (int j = columnBeginIndex; j < numOfValsToRead; j++)
                    //    //{

                    //    //}
                    //}
                }

                CP = Matrix<double>.Build.DenseOfColumnArrays(CPList.ToArray());
                //CP = Matrix<double>.Build.DenseOfColumnArrays(CPArr);                
                double[] zk = CP.QR(QRMethod.Thin).Solve(d).ToArray();
                //CP = Matrix<double>.Build.SparseOfColumnArrays(CPList.ToArray());
                //zk = Pardiso.Solve(CP, d.ToArray(), P);
                //z = Vector<double>.Build.Dense(ggg);

                zColIndex = 0;
                for (int i = 0; i < n; i++)
                {
                    if (P[i])
                    {
                        z[i] = zk[zColIndex];    //
                        zColIndex++;
                    }
                }

                //for (int i = 0; i < z.Count; i++)
                //{
                //    if (!P[i]) z[i] = 0;
                //}

                // ak je stlpec pasivny, tak ho skopcit do CPs values, aj CPs row indices, CPs row pointers uz treba zratat nanovo

                // prekopcit spravne stlpce do CPs

                //inner loop - check if any regression coefficient has turned negative
                while (z.Where((val, idx) => P[idx] == true).Any(val => val <= 0))
                {
                    innerIter++;
                    if (innerIter > iterMax)
                    {
                        throw new Exception("max pocet iteracii");
                    }

                    List<double> Q = new List<double>();

                    zColIndex = 0;

                    // get all negative regression coefficients
                    for (int i = 0; i < n; i++)
                    {
                        if (P[i] && z[i] <= 0)  // 
                        {
                            //if (z[zColIndex] <= 0) //
                            //{
                                //TODO: maybe it does not have to be a list, because we iterate through the inner loop
                                double derp = x[i] / (x[i] - z[i]);     //zColIndex
                                Q.Add(derp);
                            //}

                            //zColIndex++;
                        }
                    }

                    // find the optimal alpha to make correction with
                    double alpha = Q.Min();
                    //zColIndex = 0;

                    // find the solution that does not violate the constraints using that alpha
                    for (int i = 0; i < n; i++)
                    {
                        //if (P[i])
                        //{
                            x[i] += alpha * (z[i] - x[i]); //zColIndex
                            //zColIndex++;
                        //}
                        //else
                        //{
                        //    x[i] -= alpha * x[i];
                        //}

                        // update P and Z accordingly
                        if (Math.Abs(x[i]) < tolx && P[i])
                        {
                            A[i] = true;
                        }

                        P[i] = !A[i];
                    }

                    // recalculate for z
                    CPList.Clear();
                    //Array.Copy(Carr, CPArr, Carr.Length);

                    for (int i = 0; i < n; i++)
                    {
                        if (P[i])
                        {
                            CPList.Add(Carr[i]);
                        }

                        //if (!P[i])
                        //{
                        //    for (int j = 0; j < CPArr[i].Length; j++)
                        //    {
                        //        CPArr[i][j] = 0;
                        //    }
                        //}
                    }

                    CP = Matrix<double>.Build.DenseOfColumnArrays(CPList.ToArray());
                    //CP = Matrix<double>.Build.DenseOfColumnArrays(CPArr);
                    //CP = Matrix<double>.Build.SparseOfColumnArrays(CPList.ToArray());
                    zk = CP.QR(QRMethod.Thin).Solve(d).ToArray();
                    //zk = Pardiso.Solve(CP, d.ToArray(), P);
                    //z = Vector<double>.Build.Dense(ggg);

                    Array.Clear(z, 0, z.Length);
                    zColIndex = 0;
                    for (int i = 0; i < n; i++)
                    {
                        if (P[i])  // && z[i] <= 0
                        {
                            z[i] = zk[zColIndex];   //
                            zColIndex++;
                        }
                    }

                    //for (int i = 0; i < z.Count; i++)
                    //{
                    //    if (!P[i]) z[i] = 0;
                    //}
                }

                // calculate gradient
                //zColIndex = 0;

                //for (int i = 0; i < n; i++)
                //{
                //if (P[i])
                //{
                //x[i] = z[i];    //zColIndex
                //zColIndex++;
                //}
                //}
                Array.Copy(z, x, z.Length);
                //x = Vector<double>.Build.DenseOfVector(z);

                //resid = d - C * x;
                //w = CT * resid;

                try
                {
                    C.Multiply(x, temp.AsArray());  //TODO: may return null
                    resid = d - temp;
                    //w = CT * resid;
                    CT.Multiply(resid.AsArray(), w.AsArray());  //TODO: may return null
                }
                catch (Exception)
                {

                    throw;
                }
            }

            return Vector<double>.Build.DenseOfArray(x);
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
        [Obsolete]
        private static Vector<double> NNLS(SparseMatrix C, Vector<double> d)
        {
            //TODO: this function can be rewritten to work with sparse matrices and vectors - faster, less ram consumption, nicer

            //TODO: this can be set here as well, we dont really call that many mathnet functions
            //MathNet.Numerics.Control.UseNativeMKL();
            //MathNet.Numerics.Control.UseMultiThreading();

            //Declarations and initializations, where necessary.
            int m = C.RowCount;
            int n = C.ColumnCount;

            // related to matrix C
            SparseMatrix CT = (SparseMatrix)C.Transpose();
            double[][] Carr = new double[C.ColumnCount][];  //TODO: not effective


            for (int i = 0; i < C.ColumnCount; i++)
            {
                Carr[i] = C.Column(i);
            }

            Matrix<double> CP;
            List<double[]> CPList = new List<double[]>();
            //SparseMatrix CPs;
            //double[][] CPArr = new double[C.ColumnCount][];

            //for (int i = 0; i < C.ColumnCount; i++)
            //{
            //    CPArr[i] = new double[C.ColumnCount];
            //}

            //CP = Matrix<double>.Build.DenseOfColumnArrays(CPArr);

            // solution vectors
            Vector<double> x = Vector<double>.Build.Dense(n, 0);
            //double[] x = new double[n];
            Vector<double> temp = Vector<double>.Build.Dense(n, 0);
            Vector<double> z;

            // gradient vectors
            Vector<double> w = Vector<double>.Build.Dense(n, 0);
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

            //System.Diagnostics.Debug.WriteLine(d.Count + " " + C.RowCount + " " + C.ColumnCount + " " + x.Count);

            C.Multiply(x.AsArray(), temp.AsArray());  //TODO: may return null
            resid = d - temp;
            //w = CT * resid;
            CT.Multiply(resid.AsArray(), w.AsArray());  //TODO: may return null

            //Pardiso.Initialize(C);

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

                // pre dlzku colPointers musime najprv zratat kolko je pasivnych stlpcov
                //int pCols = 0;
                //for (int i = 0; i < n; i++)
                //{
                //    if (P[i]) pCols++;
                //}

                // kolko bude nenulovych elementov? no nie viac ako povodny pocet nenulovych

                //double[] values = new double[C.NonZerosCount];
                //int[] rowIndices = new int[C.NonZerosCount];
                //int[] colPointers = new int[pCols];

                //colPointers[0] = 0;

                for (int i = 0; i < n; i++)
                {
                    //int numOfValsToRead = C.ColumnPointers[i + 1] - C.ColumnPointers[i];
                    //int columnBeginIndex = C.ColumnPointers[i];

                    if (P[i])
                    {
                        CPList.Add(Carr[i]);                      

                        //Array.Copy(C.Values, columnBeginIndex, values, columnBeginIndex, numOfValsToRead);
                        //Array.Copy(C.RowIndices, columnBeginIndex, rowIndices, columnBeginIndex, numOfValsToRead);
                    }
                    //else    // set the rest to zero, because the columns might have been reactivated
                    //{
                    //    // TODO: zatial netreba, lebo teraz sa tie arraye robia vzdy nanovo, potom ich presuniem mimo slucku
                    //    //for (int j = columnBeginIndex; j < numOfValsToRead; j++)
                    //    //{

                    //    //}
                    //}
                }

                CP = Matrix<double>.Build.DenseOfColumnArrays(CPList.ToArray());
                //CP = Matrix<double>.Build.SparseOfColumnArrays(CPList.ToArray());
                z = CP.QR(QRMethod.Thin).Solve(d);
                //double[] ggg = Pardiso.Solve(CP, d.ToArray(), P);
                //z = Vector<double>.Build.Dense(ggg);

                //for (int i = 0; i < z.Count; i++)
                //{
                //    if (!P[i]) z[i] = 0;
                //}

                // ak je stlpec pasivny, tak ho skopcit do CPs values, aj CPs row indices, CPs row pointers uz treba zratat nanovo

                // prekopcit spravne stlpce do CPs

                //inner loop - check if any regression coefficient has turned negative
                while (z.Where((val, idx) => P[idx] == true).Any(val => val <= 0))
                {
                    innerIter++;
                    if (innerIter > iterMax)
                    {
                        throw new Exception("max pocet iteracii");
                    }

                    List<double> Q = new List<double>();

                    zColIndex = 0;

                    // get all negative regression coefficients
                    for (int i = 0; i < n; i++)
                    {
                        if (P[i])  // && z[i] <= 0
                        {
                            if (z[zColIndex] <= 0) //
                            {
                                //TODO: maybe it does not have to be a list, because we iterate through the inner loop
                                double derp = x[i] / (x[i] - z[zColIndex]);     //
                                Q.Add(derp);
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
                            x[i] += alpha * (z[zColIndex] - x[i]); //
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
                    //CP = Matrix<double>.Build.SparseOfColumnArrays(CPList.ToArray());
                    z = CP.QR(QRMethod.Thin).Solve(d);
                    //ggg = Pardiso.Solve(CP, d.ToArray(), P);
                    //z = Vector<double>.Build.Dense(ggg);

                    //for (int i = 0; i < z.Count; i++)
                    //{
                    //    if (!P[i]) z[i] = 0;
                    //}
                }

                // calculate gradient
                zColIndex = 0;

                for (int i = 0; i < n; i++)
                {
                    if (P[i])
                    {
                        x[i] = z[zColIndex];    //
                        zColIndex++;
                    }
                }
                //x = Vector<double>.Build.DenseOfVector(z);

                //resid = d - C * x;
                //w = CT * resid;

                try
                {
                    C.Multiply(x.AsArray(), temp.AsArray());  //TODO: may return null
                    resid = d - temp;
                    //w = CT * resid;
                    CT.Multiply(resid.AsArray(), w.AsArray());  //TODO: may return null
                }
                catch (Exception)
                {

                    throw;
                }
            }

            return x;
        }

        ///// <summary>
        ///// Calculate a non-negative least squares solution of C * x = d
        ///// </summary>
        ///// <remarks>
        ///// The algorithm was originally published in:
        ///// Lawson, Hanson, Solving Least Squares Problems, 1987, ISBN 978-0-89871-356-5, p. 160, Chapter 23.3
        ///// Another description of the same algorithm, but easier to understand has been published in:
        ///// Bro, De Jong, 1997, A fast non-negativity-constrained least squares algorithm, J. Chemom. 11, pp. 393-401
        ///// </remarks>
        ///// <param name="C">Matrix describing the model.</param>
        ///// <param name="d">Vector with observation values.</param>
        ///// <returns>MathNet vector with the solution.</returns>
        //private static Vector<double> NNLS(SparseMatrix C, Vector<double> d)
        //{
        //    //TODO: this needs to be set at the start of all calculations, right after Isotopefitter is called
        //    //MathNet.Numerics.Control.UseNativeMKL();
        //    //MathNet.Numerics.Control.UseMultiThreading();

        //    //Declarations and initializations, where necessary.
        //    int m = C.RowCount;
        //    int n = C.ColumnCount;

        //    // related to matrix C
        //    //Matrix<double> CT = C.Transpose();
        //    SparseMatrix CT = (SparseMatrix)C.Transpose();

        //    //double[][] Carr = C.ToColumnArrays();
        //    //TODO: here
        //    //double[][] Carr = C.

        //    Matrix<double> CP;
        //    List<double[]> CPList = new List<double[]>();

        //    // solution vectors
        //    Vector<double> x = Vector<double>.Build.Dense(n, 0);
        //    Vector<double> z;

        //    // gradient vectors
        //    Vector<double> w;
        //    Vector<double> wz = Vector<double>.Build.Dense(n, 0);

        //    // active and passive sets
        //    bool[] P = new bool[n];
        //    bool[] A = new bool[n];

        //    for (int i = 0; i < n; i++)
        //    {
        //        P[i] = false;
        //        A[i] = true;
        //    }

        //    // helper variables
        //    Vector<double> resid;
        //    int outerIter = 0;
        //    int innerIter = 0;
        //    int iterMax = 3 * n;
        //    int zColIndex = 0;
        //    double tolx = 10 * MathNet.Numerics.Precision.DoublePrecision * n * C.L1Norm(); // that is how octave does it

        //    // calculation starts here

        //    System.Diagnostics.Debug.WriteLine(d.Count + " " + C.RowCount + " " + C.ColumnCount + " " + x.Count);

        //    resid = d - C * x;
        //    w = CT * resid;

        //    /* 
        //     * Main loop.
        //     * First condition checks whether there are still design parameters in the active set.
        //     * The second condition first filters out those elements from vector w, that correspond to active design parameters.
        //     * Then, those filtered values are compared to convergence tolerance.
        //     * In human language: In the end, the active set should be either empty, or the corresponding w elements must be close to zero.
        //     */
        //    while (A.Any(b => b == true) && w.Where((val, idx) => A[idx] == true).Any(val => val > tolx))
        //    {
        //        outerIter += 1;


        //        //TODO: check if wz initialization can be implemented using LINQ, something like wz.Where((val, idx) => P[idx] == false).Select(val => val = double.MinValue);
        //        for (int i = 0; i < n; i++)
        //        {
        //            System.Diagnostics.Debug.Assert(P[i] != A[i]);

        //            if (P[i])
        //            {
        //                wz[i] = double.MinValue;
        //            }
        //            else    // that means Z[i] is true
        //            {
        //                wz[i] = w[i];
        //            }
        //        }

        //        // we move the design parameter with the highest gradient (w) to the passive set
        //        int t = wz.MaximumIndex();
        //        P[t] = true;
        //        A[t] = false;

        //        // find intermediary solution (z)
        //        CPList.Clear();

        //        for (int i = 0; i < n; i++)
        //        {
        //            if (P[i])
        //            {
        //                CPList.Add(Carr[i]);
        //            }
        //        }

        //        CP = Matrix<double>.Build.DenseOfColumnArrays(CPList.ToArray());
        //        z = CP.QR(QRMethod.Thin).Solve(d);

        //        //inner loop - check if any regression coefficient has turned negative
        //        while (z.Where((val, idx) => P[idx] == true).Any(val => val <= 0))
        //        {
        //            innerIter++;
        //            if (innerIter > iterMax)
        //            {
        //                Console.WriteLine("max pocet iteracii");
        //            }

        //            List<double> Q = new List<double>();

        //            zColIndex = 0;

        //            // get all negative regression coefficients
        //            for (int i = 0; i < n; i++)
        //            {
        //                if (P[i])
        //                {
        //                    if (z[zColIndex] <= 0)
        //                    {
        //                        //TODO: maybe it does not have to be a list, because we iterate through the inner loop
        //                        Q.Add(x[i] / (x[i] - z[zColIndex]));
        //                    }

        //                    zColIndex++;
        //                }
        //            }

        //            // find the optimal alpha to make correction with
        //            double alpha = Q.Min();
        //            zColIndex = 0;

        //            // find the solution that does not violate the constraints using that alpha
        //            for (int i = 0; i < n; i++)
        //            {
        //                if (P[i])
        //                {
        //                    x[i] += alpha * (z[zColIndex] - x[i]);
        //                    zColIndex++;
        //                }
        //                else
        //                {
        //                    x[i] -= alpha * x[i];
        //                }

        //                // update P and Z accordingly
        //                if (Math.Abs(x[i]) < tolx && P[i])
        //                {
        //                    A[i] = true;
        //                }

        //                P[i] = !A[i];
        //            }

        //            // recalculate for z
        //            CPList.Clear();

        //            for (int i = 0; i < n; i++)
        //            {
        //                if (P[i])
        //                {
        //                    CPList.Add(Carr[i]);
        //                }
        //            }

        //            CP = Matrix<double>.Build.DenseOfColumnArrays(CPList.ToArray());
        //            z = CP.QR(QRMethod.Thin).Solve(d);
        //        }

        //        // calculate gradient
        //        zColIndex = 0;

        //        for (int i = 0; i < n; i++)
        //        {
        //            if (P[i])
        //            {
        //                x[i] = z[zColIndex];
        //                zColIndex++;
        //            }
        //        }

        //        resid = d - C * x;
        //        w = CT * resid;
        //    }

        //    return x;
        //}
    }
}
