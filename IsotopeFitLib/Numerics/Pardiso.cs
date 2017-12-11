using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Security;
using System.Runtime.InteropServices;

using CSparse.Double;

namespace IsotopeFit
{
    public sealed class Pardiso
    {
        static IntPtr[] pt;
        static int maxfct = 1;
        static int mnum = 1;
        static int mtype = 11; // real nonsymmetric
        static int phase = 0;
        static int n = 0;
        static int nrhs = 1;
        static int[] iparm;
        static int msglvl = 1;
        static int err = 0;

        static double[] ddum; /* Double dummy */
        static int[] idum; /* Integer dummy. */

        private Pardiso() { }

        public static void Initialize(SparseMatrix A)
        {
            // matrix comes in the compressed column format, pardiso wants compressed row.
            SparseMatrix At = (SparseMatrix)A.Transpose();

            // then, the internal arrays contain:
            double[] values = At.Values;
            int[] columnIndices = At.RowIndices;
            int[] rowPointers = At.ColumnPointers;

            pt = new IntPtr[64];
            phase = 12;
            n = A.RowCount;

            InitIparm();

            idum = new int[n];
            ddum = new double[n];

            //CallSolver(pt, ref maxfct, ref mnum, ref mtype, ref phase, ref n, values, rowPointers, columnIndices, idum, ref nrhs, iparm, ref msglvl, ddum, ddum, ref err);

            if (err != 0) throw new Exception("Pardiso initialization failed.");
        }

        /// <summary>
        /// Solves the linear equation system using Intel MKL Pardiso solver.
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="P"></param>
        public static double[] Solve(MathNet.Numerics.LinearAlgebra.Matrix<double> A, double[] b, bool[] P)
        {
            // matrix comes in the compressed column format, pardiso wants compressed row.
            MathNet.Numerics.LinearAlgebra.Double.SparseMatrix At = (MathNet.Numerics.LinearAlgebra.Double.SparseMatrix)A.Transpose();

            // then, the internal arrays contain:
            double[] values = (A.Storage as MathNet.Numerics.LinearAlgebra.Storage.SparseCompressedRowMatrixStorage<double>).Values;
            int[] columnIndices = (A.Storage as MathNet.Numerics.LinearAlgebra.Storage.SparseCompressedRowMatrixStorage<double>).ColumnIndices;
            int[] rowPointers = (A.Storage as MathNet.Numerics.LinearAlgebra.Storage.SparseCompressedRowMatrixStorage<double>).RowPointers;

            //IntPtr[] pt = new IntPtr[64];
            //int n = 0;  //A.RowCount; // TODO: number of equations (non-zero rows?)

            int[] perm = new int[P.Length];

            for (int i = 0; i < P.Length; i++)
            {
                if (P[i])
                {
                    //perm[i] = 1;
                    //n++;
                }
                else
                {
                    //b[i] = 0;
                }
            }


            //n = A.RowCount;

            //int nrhs = 1;

            //int[] iparm = new int[64];
            //for (int i = 0; i < 64; i++)
            //{
            //    iparm[i] = 0;
            //}
            //iparm[0] = 1; /* No solver default */
            //iparm[1] = 0; // minimum degree
            //              /* Numbers of processors, value of OMP_NUM_THREADS */
            //iparm[2] = 0; // why?
            //iparm[3] = 0; /* No iterative-direct algorithm */
            //iparm[4] = 0; // no permutation input
            //iparm[5] = 0; /* Write solution into x */
            //iparm[6] = 0; // output, number of refinement steps
            //iparm[7] = 0; /* Max numbers of iterative refinement steps */
            //iparm[8] = 0; // reserved
            //iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
            //iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
            //iparm[11] = 0; // here we can set if the A matrix gets intenally transposed?
            //iparm[12] = 0; // Maximum weighted matching algorithm is switched-off. TODO: switch to 1 if accuracy is low
            //iparm[13] = 0; /* Output: Number of perturbed pivots */

            //iparm[17] = 0; /* Output: Number of nonzeros in the factor LU */
            //iparm[18] = 0; /* Output: Mflops for LU factorization */
            //iparm[19] = 0; /* Output: Numbers of CG Iterations */

            //iparm[23] = 0; // parallel setting for factorization (maybe 1)
            //iparm[24] = 0; // parallel solve
            //iparm[26] = 1; // indices check

            //iparm[30] = 3; // take advantage of the right-hand side sparsity
            //iparm[33] = 0; // CNR mode
            //iparm[34] = 1; // zero based indexing
            //iparm[35] = 0;

            //iparm[59] = 0; // in core

            //int maxfct = 1;
            //int mnum = 1;
            //int mtype = 11; // real nonsymmetric
            //int phase = 13; // Analysis, numerical factorization, solve, iterative refinement
            //int msglvl = 0;
            //double[] sol = new double[b.Length];
            //int err = 0;

            //CallSolver(pt, ref maxfct, ref mnum, ref mtype, ref phase, ref n, values, rowPointers, columnIndices, idum, ref nrhs, iparm, ref msglvl, ddum, ddum, ref err);

            //iparm[4] = 0;
            //iparm[30] = 0;
            //iparm[33] = 0;
            //phase = 12;

            //CallSolver(pt, ref maxfct, ref mnum, ref mtype, ref phase, ref n, values, rowPointers, columnIndices, idum, ref nrhs, iparm, ref msglvl, ddum, ddum, ref err);

            //InitIparm();

            pt = new IntPtr[64];
            //phase = 12;
            n = A.RowCount;

            InitIparm();

            //iparm[30] = 3;
            phase = 13;

            double[] sol = new double[b.Length];


            CallSolver(pt, ref maxfct, ref mnum, ref mtype, ref phase, ref n, values, rowPointers, columnIndices, perm, ref nrhs, iparm, ref msglvl, b, sol, ref err);

            return sol;
        }

        /// <summary>
        /// Intel MKL pardiso wrapper method. Calculates the solution of a set of sparse linear equations A*X = B with single or multiple right-hand sides.
        /// </summary>
        /// <param name="handle">Handle to internal data structure. Also on output.</param>
        /// <param name="maxfct">Maximum number of factors to be stored (default 1).</param>
        /// <param name="mnum">Which matrix to factorize (default 1).</param>
        /// <param name="mtype">Defines the matrix type, influences pivoting. (see MKL reference for values)</param>
        /// <param name="phase">Execution control. (see MKL reference for values)</param>
        /// <param name="n">Number of equations in the system.</param>
        /// <param name="a">Values of the coefficient matrix A, corresponding to the indices in ja. (CSR format)</param>
        /// <param name="ia">Row pointers of A matrix values. (CSR format)</param>
        /// <param name="ja">Column indices of A matrix values. (CSR format)</param>
        /// <param name="perm">Holds the permutation vector. Depends on iparm settings (see MKL reference). Also on output.</param>
        /// <param name="nrhs">Number of right hand sides.</param>
        /// <param name="iparm">Array[64] to pass settings and retrieve information. Also on output.</param>
        /// <param name="msglvl">Message output level. (0 - no messages)</param>
        /// <param name="b">Array containing the right hand side. On output contains the solution if iparm is set correspondingly.</param>
        /// <param name="x">Contains solution vector on output.</param>
        /// <param name="error">Error code. (see MKL reference for values)</param>
        /// <returns>Dunno what this is. MKL reference says it is void.</returns>
        internal static int CallSolver(IntPtr[] handle,
            ref int maxfct, ref int mnum,
            ref int mtype, ref int phase, ref int n,
            double[] a, int[] ia, int[] ja, int[] perm,
            ref int nrhs, int[] iparm, ref int msglvl,
            double[] b, double[] x, ref int error)
        {
            return PardisoNative.pardiso(handle,
                ref maxfct, ref mnum, ref mtype, ref phase, ref n,
                a, ia, ja, perm, ref nrhs, iparm, ref msglvl,
                b, x, ref error);
        }

        private static void InitIparm()
        {
            iparm = new int[64];
            for (int i = 0; i < 64; i++)
            {
                iparm[i] = 0;
            }
            //iparm[0] = 1; /* No solver default */
            //iparm[1] = 2; // nested dissection
            //              /* Numbers of processors, value of OMP_NUM_THREADS */
            //iparm[2] = 0; // why?
            //iparm[3] = 0; /* No iterative-direct algorithm */
            //iparm[4] = 0; // no permutation input
            //iparm[5] = 0; /* Write solution into x */
            //iparm[6] = 0; // output, number of refinement steps
            //iparm[7] = 0; /* Max numbers of iterative refinement steps */
            //iparm[8] = 0; // reserved
            //iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
            //iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
            //iparm[11] = 0; // here we can set if the A matrix gets intenally transposed?
            //iparm[12] = 1; // Maximum weighted matching algorithm is switched-off. TODO: switch to 1 if accuracy is low
            //iparm[13] = 0; /* Output: Number of perturbed pivots */

            //iparm[17] = 0; /* Output: Number of nonzeros in the factor LU */
            //iparm[18] = 0; /* Output: Mflops for LU factorization */
            //iparm[19] = 0; /* Output: Numbers of CG Iterations */

            //iparm[23] = 0; // parallel setting for factorization (maybe 1)
            //iparm[24] = 0; // parallel solve
            //iparm[26] = 1; // indices check

            //iparm[30] = 0; // take advantage of the right-hand side sparsity
            //iparm[33] = 0; // CNR mode
            iparm[34] = 1; // zero based indexing
            //iparm[35] = 0;

            //iparm[59] = 0; // in core
        }
    }

    /** Pardiso native declarations */
    [SuppressUnmanagedCodeSecurity]
    internal sealed class PardisoNative
    {
        private PardisoNative() { }

        [DllImport("mkl.dll", CallingConvention = CallingConvention.Cdecl,
             ExactSpelling = true, SetLastError = false)]
        internal static extern int pardiso(
            [In, Out] IntPtr[] handle,
            ref int maxfct,
            ref int mnum,
            ref int mtype, ref int phase, ref int n,
            [In] double[] a, [In] int[] ia, [In] int[] ja, [In] int[] perm,
            ref int nrhs, [In, Out] int[] iparm, ref int msglvl,
            [In, Out] double[] b, [Out] double[] x, ref int error);
    }
}
