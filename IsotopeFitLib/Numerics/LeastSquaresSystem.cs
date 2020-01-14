using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.CompilerServices;

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
        public LeastSquaresSystem(SparseMatrix desMat, MathNet.Numerics.LinearAlgebra.Double.SparseVector obsVec)
        {
            DesignMatrix = desMat;
            ObservationVector = obsVec;
        }

        internal SparseMatrix DesignMatrix { get; private set; }
        internal SparseMatrix DesignMatrixR { get; set; }
        //internal bool[] DMFitmask { get; set; }
        internal MathNet.Numerics.LinearAlgebra.Double.SparseVector ObservationVector { get; private set; }
        internal Vector<double> ObservationVectorR { get; set; }
        internal int[] ColumnOrdering { get; set; }
        public double[] Solution { get; private set; }
        //internal Vector<double> FittedSpectrum { get; private set; }
        public double[] SolutionError { get; private set; }

        /// <summary>
        /// Method that solves the system.
        /// </summary>
        /// <remarks>
        /// Partitions the problem by non-overlapping diagonal elements and non-overlapping blocks.
        /// </remarks>
        public async Task Solve2Async()
        {
            double[] values = DesignMatrixR.Values;
            int[] rowIndices = DesignMatrixR.RowIndices;
            int[] columnPointers = DesignMatrixR.ColumnPointers;

            SparseMatrix AT = (SparseMatrix)DesignMatrixR.Transpose();

            double[] valuesT = AT.Values;
            int[] columnIndicesT = AT.RowIndices;
            int[] rowPointersT = AT.ColumnPointers;

            double[] colCounts = new double[columnPointers.Length - 1];
            double[] rowCounts = new double[rowPointersT.Length - 1];

            List<int> cutCoordinates = new List<int>();

            for (int i = 0; i < colCounts.Length; i++)
            {
                colCounts[i] = columnPointers[i + 1] - columnPointers[i];
                rowCounts[i] = rowPointersT[i + 1] - rowPointersT[i];
            }

            // it is divided into two loops, because we nee the i+1st element of the array
            for (int i = 0; i < colCounts.Length; i++)
            {
                colCounts[i] = columnPointers[i + 1] - columnPointers[i];
                rowCounts[i] = rowPointersT[i + 1] - rowPointersT[i];

                if (rowCounts[i] == 1 && colCounts[i] == 1)     // partitioning 1
                {
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

            // calculate the block sizes
            int[] blSizes = new int[cutCoordinates.Count];
            blSizes[0] = cutCoordinates[0] + 1;

            for (int i = 1; i < blSizes.Length; i++)
            {
                blSizes[i] = cutCoordinates[i] - cutCoordinates[i - 1]; //TODO: currently it is block sizes and blocks ends - maybe block beginnins would be nicer?
            }

            // calculate the non-overlaping cluster abundances
            double[] abd = new double[DesignMatrixR.ColumnCount];
            double[] abdErrors = new double[DesignMatrixR.ColumnCount];

            List<Task> nnlsTaskList = new List<Task>();
            List<Task> errorCalcTaskList = new List<Task>();

            for (int i = 0; i < cutCoordinates.Count; i++)
            {
                if (blSizes[i] == 1)
                {
                    abd[cutCoordinates[i]] = ObservationVectorR[cutCoordinates[i]] / DesignMatrixR.At(cutCoordinates[i], cutCoordinates[i]);

                    // TODO: error estimation
                    // TODO: take the cutCoordinates[i]-th column of the original design matrix and multiply by the abundance. then proceed according to the formula
                    // TODO: Wrong! We need to take into account the column ordering as well!!!
                    abdErrors[cutCoordinates[i]] = await CalculateError(DesignMatrix, cutCoordinates[i], abd[cutCoordinates[i]]);  
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
                    Array.Copy(ObservationVectorR.ToArray(), start, obs, 0, end - start);

                    //Console.WriteLine("starting task {0}", i);
                    nnlsTaskList.Add(Task.Run(() =>
                    {
                        double[] sol = NNLS3(C, obs);
                        Array.Copy(sol, 0, abd, start, sol.Length);

                        //errorCalcTaskList.Add(Task.Run(() => CalculateError2(C, start, sol)));
                        //CalculateError2(C, start, sol);
                    }));
                    // error estimation

                    // TODO: temporary, to ease debugging
                    //nnlsTaskList.Last().Wait();
                    //errorCalcTaskList.Last().Wait();
                }
            }

            Task.WaitAll(nnlsTaskList.ToArray());

            Solution = abd;            

            Solution = Enumerable.Zip(ColumnOrdering, Solution, (idx, val) => new { idx, val }).OrderBy(v => v.idx).Select(v => v.val).ToArray();
            //SolutionError = Enumerable.Zip(ColumnOrdering, SolutionError, (idx, val) => new { idx, val }).OrderBy(v => v.idx).Select(v => v.val).ToArray();

            //Task.WaitAll(errorCalcTaskList.ToArray());
            //CalculateError3(DesignMatrix, DesignMatrixR, ObservationVector, Solution);



        }

        //TODO: this is the primordial form of the error calculation function. It will allow for calculation of fit errors on demand from the GUI.
        //TODO: at the moment, this is only for performance testing purposes
        internal void CalcErr()
        {
            CalculateError3(DesignMatrix, DesignMatrixR, ObservationVector, Solution);
        }

        private Task<double> CalculateError(SparseMatrix designMatrix, int idx, double abundance)
        {
            return Task.Run((Func<double>)(() =>
            {
                // remap according to column reorderings
                int orderIdx = ColumnOrdering.ToList().FindIndex((Predicate<int>)(a => a == idx));

                double[] calcSignal = new double[designMatrix.ColumnPointers[orderIdx + 1] - designMatrix.ColumnPointers[orderIdx]];
                int j = 0;
                double diffSqSum = 0;

                for (int i = designMatrix.ColumnPointers[orderIdx]; i < designMatrix.ColumnPointers[orderIdx + 1]; i++)
                {
                    calcSignal[j] = designMatrix.Values[i] * abundance;
                    diffSqSum += Math.Pow(calcSignal[j] - ObservationVector.At(designMatrix.RowIndices[i]), 2);
                    j++;
                }
                
                double sSqInv = 1d / Math.Pow(DesignMatrixR.Values[DesignMatrixR.ColumnPointers[idx]], 2);

                return 1.96d * Math.Sqrt(sSqInv * diffSqSum / (calcSignal.Length - 1));
            }));
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
            Solution = NNLS3(DesignMatrix, ObservationVector.ToArray());

            // reorder the solution according to the ColumnOrdering information from sparse QR factorization
            Solution = Enumerable.Zip(ColumnOrdering, Solution, (idx, val) => new { idx, val }).OrderBy(v => v.idx).Select(v => v.val).ToArray();
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
        private static double[] NNLS3(SparseMatrix C, double[] d)
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
            SparseMatrix CP;

            // solution vectors
            double[] x = new double[n];
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
            Vector<double> dv = Vector<double>.Build.DenseOfArray(d);
            Vector<double> resid;
            int outerIter = 0;
            int innerIter = 0;
            int iterMax = 3 * n;
            int zColIndex = 0;
            double tolx = 10 * MathNet.Numerics.Precision.DoublePrecision * n * C.L1Norm(); // that is how octave does it

            // calculation starts here
            C.Multiply(x, temp.AsArray());  //TODO: may return null
            resid = dv - temp;
            //w = CT * resid;
            CT.Multiply(resid.AsArray(), w.AsArray());  //TODO: may return null

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
                        if (P[i] && z[i] <= 0)
                        {
                            //TODO: maybe it does not have to be a list, because we iterate through the inner loop
                            double derp = x[i] / (x[i] - z[i]);
                            Q.Add(derp);
                        }
                    }

                    // find the optimal alpha to make correction with
                    double alpha = Q.Min();

                    // find the solution that does not violate the constraints using that alpha
                    for (int i = 0; i < n; i++)
                    {
                        x[i] += alpha * (z[i] - x[i]);

                        // update P and A accordingly
                        if (Math.Abs(x[i]) < tolx && P[i])
                        {
                            A[i] = true;
                        }

                        P[i] = !A[i];
                    }
                    
                    CP = BuildCP(C, P);
                    zk = SPQR.Solve(CP, d.ToArray());

                    Array.Clear(z, 0, z.Length);
                    zColIndex = 0;
                    for (int i = 0; i < n; i++)
                    {
                        if (P[i])
                        {
                            z[i] = zk[zColIndex];
                            zColIndex++;
                        }
                    }
                }

                // calculate gradient
                Array.Copy(z, x, z.Length);

                C.Multiply(x, temp.AsArray());  //TODO: may return null
                resid = dv - temp;
                //w = CT * resid;
                CT.Multiply(resid.AsArray(), w.AsArray());  //TODO: may return null
            }

            return x;
        }

        /// <summary>
        /// Creates a new matrix from selected columns of original matrix.
        /// </summary>
        /// <param name="C">Original matrix.</param>
        /// <param name="P">Bool array designating which columns should comprise the new matrix.</param>
        /// <returns>New matrix built from the selected columns.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static SparseMatrix BuildCP(SparseMatrix C, bool[] P)
        {
            // count how many active columns we have
            int pCols = 0;
            int nnzCount = 0;
            int columnNnzCount = 0;

            for (int i = 0; i < C.RowCount; i++)
            {
                if (P[i])
                {
                    pCols++;
                    nnzCount += C.ColumnPointers[i + 1] - C.ColumnPointers[i];
                }
            }

            // construct the matrix arrays
            double[] Avalues = new double[nnzCount];
            int[] ArowIndices = new int[nnzCount];
            int[] AcolPointers = new int[pCols + 1];

            pCols = 0;
            for (int i = 0; i < C.RowCount; i++)
            {
                if (P[i])
                {
                    columnNnzCount = C.ColumnPointers[i + 1] - C.ColumnPointers[i];

                    AcolPointers[pCols + 1] = AcolPointers[pCols] + columnNnzCount;                    

                    Array.Copy(C.Values, C.ColumnPointers[i], Avalues, AcolPointers[pCols], columnNnzCount);
                    Array.Copy(C.RowIndices, C.ColumnPointers[i], ArowIndices, AcolPointers[pCols], columnNnzCount);

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
    }
}
