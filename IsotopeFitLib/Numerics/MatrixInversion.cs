using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using CSparse.Double;
using MathNet.Numerics.LinearAlgebra;
using System.Diagnostics;

namespace IsotopeFit.Numerics
{
    public static class MatrixInversion
    {
        static SparseMatrix A; // Original Matrix
        static SparseMatrix ATransposed; // Transposed Original Matrix
        static int n;
        static double tolx; 

        static Col[] ColArray;
        static double[] invVal;
        static int[] invRwIdx;
        static int[] invClnIdx;
        static SparseMatrix invA;
        
        //static Stopwatch sw = new Stopwatch();

        public static SparseMatrix Inverse(SparseMatrix a)
        {
            A = a;
            n = A.RowCount;
            ColArray = new Col[n];
            tolx = 10 * MathNet.Numerics.Precision.DoublePrecision * (n - 1) * A.L1Norm(); // Tolerance

            //sw.Start();

            Parallel.For(0, n, InvertColumn); // columns of inversed matrix do not influence each other

            int nonZeros = 0; // total number of non-zero values in inversed matrix
            for (int i = 0; i < ColArray.Length; i++)
            {
                nonZeros += ColArray[i].Val.Count;
            }

            invVal = new double[nonZeros];
            invRwIdx = new int[nonZeros];
            invClnIdx = new int[n + 1];
            int cc = 0;

            for (int i = 0; i < n; i++)
            {
                ColArray[i].Val.CopyTo(invVal, cc);
                ColArray[i].RwIdx.CopyTo(invRwIdx, cc);
                cc += ColArray[i].Val.Count;
                invClnIdx[i + 1] = cc;
            }

            invA = new SparseMatrix(n, n)
            {
                Values = invVal,
                RowIndices = invRwIdx,
                ColumnPointers = invClnIdx
            };

            return invA;

            //sw.Stop();
            //Console.WriteLine("time (ms): " + sw.ElapsedMilliseconds);
            //Console.ReadKey();

            //for (int i = 0; i <= n; i++)
            //{
            //    for (int j = 0; j <= n; j++)
            //    {
            //        Console.Write(A.At(i, j) + "\t");
            //    }
            //    Console.WriteLine();
            ////}
            //Console.WriteLine("\n\nInv(A):");
            //for (int i = 0; i < n; i++)
            //{
            //    for (int j = 0; j < n; j++)
            //    {
            //        Console.Write(invA.At(i, j) + "\t");
            //    }
            //    Console.WriteLine();
            //}
            //Console.ReadKey();
        }

        static void InvertColumn(int i, ParallelLoopState pls)
        {
            int colSize = n - i - 1;

            double v;
            int idx;
            int rw;
            int nonZeros;

            double Xval;
            double subSum;
            double[] XVectro = new double[n];

            ColArray[colSize] = new Col();

            ColArray[colSize].Val.Add(1 / A.At(colSize, colSize));
            ColArray[colSize].RwIdx.Add(colSize);
            XVectro[colSize] = ColArray[colSize].Val[0];

            for (int j = 0; j <= colSize; j++)
            {
                idx = ATransposed.ColumnPointers[colSize - j];
                nonZeros = ATransposed.ColumnPointers[colSize + 1 - j] - idx;

                subSum = 0;

                for (int k = 1; k < nonZeros; k++)
                {
                    v = ATransposed.Values[idx + k];
                    rw = ATransposed.RowIndices[idx + k];

                    Xval = XVectro[rw];
                    if (Math.Abs(Xval) > tolx)
                        subSum += v * Xval;
                }
                if (subSum != 0)
                {
                    ColArray[colSize].Val.Add((-1) * subSum / ATransposed.Values[idx]);
                    ColArray[colSize].RwIdx.Add(ATransposed.RowIndices[idx]);
                    XVectro[ATransposed.RowIndices[idx]] = (-1) * subSum / ATransposed.Values[idx];
                }
            }
            ColArray[colSize].Val.Reverse();
            ColArray[colSize].RwIdx.Reverse();
        }
    }

    internal class Col
    {
        internal Col()
        {
            Val = new List<double>();
            RwIdx = new List<int>();
        }

        internal List<double> Val;
        internal List<int> RwIdx;
    }
}
