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

        public static SparseMatrix Calculate(SparseMatrix A)
        {
            IntPtr[] handles = new IntPtr[3];
            int ordering = 5;   // 5 - AMD mimimum degree, 6 - metis mindeg, 
            double tolerance = 1e-9;

            int rows = A.RowCount;
            int cols = A.ColumnCount;
            int nzCount = A.NonZerosCount;
            double[] values = A.Values;
            long[] rowIndices = Array.ConvertAll<int, long>(A.RowIndices, a => a);
            long[] colPointers = Array.ConvertAll<int, long>(A.ColumnPointers, a => a);

            int handle = SparseQR(handles,
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
    }
}
