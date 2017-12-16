using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

using CSparse.Double;

namespace IsotopeFit
{
    public sealed class SPQR
    {
        private SPQR() { }

        [DllImport("SPQR.dll")]
        private static extern void SparseQRDispose(IntPtr[] handle);

        [DllImport("SPQR.dll")]
        private static extern int SparseQR(
            IntPtr[] handle,
            [In] int ord, [In] double tol,
            [In] int rows, [In] int cols, [In] int nzCount, [In] double[] vals, [In] long[] rowInd, [In] long[] colPtr,
            out int Rrows, out int Rcols, out int RnzCount, out IntPtr Rvals, out IntPtr RrowInd, out IntPtr RcolPtr, out IntPtr ROrdering);

        [DllImport("SPQR.dll")]
        private static extern void SparseSolveDispose(IntPtr[] handle);

        [DllImport("SPQR.dll")]
        private static extern int SparseSolve(
            IntPtr[] handle,
            [In] int rows, [In] int cols, [In] int nzCount, [In] double[] vals, [In] long[] rowInd, [In] long[] colPtr,
            double[] bVals,
            out IntPtr xVals);

        /// <summary>
        /// Calculates the QR factorization of a sparse matrix.
        /// </summary>
        /// <remarks>
        /// Discards the calculated values of the Q matrix, in order to be faster and more memory efficient.
        /// </remarks>
        /// <param name="A">Sparse matrix to be factorized.</param>
        /// <returns>Upper triangular factor R, in compressed sparse column format.</returns>
        public static SparseMatrix QR(SparseMatrix A)
        {
            IntPtr[] handles = new IntPtr[3];
            int ordering = 5;
            double tolerance = 1e-9;

            int rows = A.RowCount;
            int cols = A.ColumnCount;
            int nzCount = A.NonZerosCount;
            double[] values = A.Values;
            long[] rowIndices = Array.ConvertAll<int, long>(A.RowIndices, a => a);
            long[] colPointers = Array.ConvertAll<int, long>(A.ColumnPointers, a => a);

            int status = SparseQR(handles,
                ordering, tolerance,
                rows, cols, nzCount, values, rowIndices, colPointers,
                out int Rrows, out int Rcols, out int RnzCount, out IntPtr Rvals, out IntPtr RrowInd, out IntPtr RcolPtr, out IntPtr RorderingPtr);

            double[] RvaluesArr = new double[RnzCount];
            Int64[] RrowIndArr = new Int64[RnzCount];
            Int64[] RcolPtrArr = new Int64[Rcols + 1];
            Int64[] RorderingArr = new Int64[Rcols];

            Marshal.Copy(Rvals, RvaluesArr, 0, RnzCount);
            Marshal.Copy(RrowInd, RrowIndArr, 0, RnzCount);
            Marshal.Copy(RcolPtr, RcolPtrArr, 0, Rcols + 1);
            Marshal.Copy(RorderingPtr, RorderingArr, 0, Rcols);

            SparseQRDispose(handles);

            SparseMatrix R = new SparseMatrix(Rrows, Rcols)
            {
                Values = RvaluesArr,
                RowIndices = Array.ConvertAll(RrowIndArr, a => (int)a),
                ColumnPointers = Array.ConvertAll(RcolPtrArr, a => (int)a),
                ColumnOrdering = Array.ConvertAll(RorderingArr, a => (int)a)
            };

            return R;
        }

        /// <summary>
        /// Solves the sparse linear equation system A * x = b.
        /// </summary>
        /// <param name="A">Sparse matrix of coefficients.</param>
        /// <param name="b">Array of observation values.</param>
        /// <returns>Array containing the solution of the linear equation system.</returns>
        public static double[] Solve(SparseMatrix A, double[] b)
        {
            IntPtr[] handles = new IntPtr[4];

            int rows = A.RowCount;
            int cols = A.ColumnCount;
            int nzCount = A.NonZerosCount;
            double[] values = A.Values;
            long[] rowIndices = Array.ConvertAll<int, long>(A.RowIndices, a => a);
            long[] colPointers = Array.ConvertAll<int, long>(A.ColumnPointers, a => a);

            int status = SparseSolve(handles, rows, cols, nzCount, values, rowIndices, colPointers, b, out IntPtr xVals);

            double[] x = new double[rows];
            Marshal.Copy(xVals, x, 0, status);

            SparseSolveDispose(handles);

            return x;
        }
    }
}
