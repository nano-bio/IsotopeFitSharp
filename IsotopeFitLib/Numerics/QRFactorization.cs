using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;

namespace IsotopeFit.Numerics
{
    public static partial class Algorithm
    {
        /// <summary>
        /// Calculates QR factorization of a least squares system.
        /// </summary>
        /// <param name="lss">Least squares system object to be QR factorized.</param>
        /// <returns>QR factorized least squares system object.</returns>
        public static LeastSquaresSystem LeaSqrSparseQR(LeastSquaresSystem lss)
        {
            Matrix<double> M = Matrix<double>.Build.SparseOfMatrix(lss.SystemMatrix);
            Vector<double> v = Vector<double>.Build.DenseOfVector(lss.ObservationVector);

            Givens gp;
            double x, y;

            for (int i = 0; i < M.ColumnCount; i++)
            {
                for (int j = M.RowCount - 1; j > i; j--)
                {
                    if (M.At(j, i) == 0) continue;

                    gp = GivensParams(M.At(i, i), M.At(j, i));

                    for (int k = i; k < M.ColumnCount; k++)
                    {
                        x = M.At(i, k);
                        y = M.At(j, k);

                        M.At(i, k, gp.c * x + gp.s * y);
                        M.At(j, k, gp.c * y - gp.s * x);

                        //TODO: add zero-coercing
                    }

                    x = v.At(i);
                    y = v.At(j);

                    //TODO: this must be tested, if it is done correctly
                    v.At(i, gp.c * x + gp.s * y);
                    v.At(j, gp.c * y - gp.s * x);
                }
            }

            //TODO: this cutting, especially for the vector should maybe be checked for non-zero values before the cutting line
            return new LeastSquaresSystem(M.SubMatrix(0, M.ColumnCount, 0, M.ColumnCount), v.SubVector(0, M.ColumnCount));
        }

        /// <summary>
        /// Calculates QR factorization of the sparse matrix M.
        /// </summary>
        /// <remarks>
        /// Method calculates the upper triangular matrix R without explicitly forming the Q matrix.
        /// It is meant primarily for sparse matrices.
        /// </remarks>
        /// <param name="M">Matrix to be QR factorized.</param>
        public static Matrix<double> SparseQR(Matrix<double> M)
        {
            Givens gp;
            double x, y;

            for (int i = 0; i < M.ColumnCount; i++)
            {
                for (int j = M.RowCount - 1; j > i; j--)
                {
                    if (M.At(j, i) == 0) continue;

                    gp = GivensParams(M.At(i, i), M.At(j, i));

                    for (int k = i; k < M.ColumnCount; k++)
                    {
                        x = M.At(i, k);
                        y = M.At(j, k);

                        M.At(i, k, gp.c * x + gp.s * y);
                        M.At(j, k, -gp.s * x + gp.c * y);

                        //TODO: add zero-coercing
                    }
                }
            }

            return M.SubMatrix(0, M.ColumnCount, 0, M.ColumnCount);
        }

        /// <summary>
        /// Calculates Givens rotation elements.
        /// </summary>
        /// <param name="v1">Matrix element to be rotatet into.</param>
        /// <param name="v2">Matrix element to be zeroed.</param>
        /// <returns>Struct containing the Givens rotation elements.</returns>
        private static Givens GivensParams(double v1, double v2)
        {
            Givens gp = new Givens();

            if (v2 == 0)
            {
                gp.c = 1;
                gp.s = 0;
            }
            else if (v1 == 0)
            {
                gp.c = 0;
                gp.s = 1;
            }
            else
            {
                //TODO: test if putting these operations is more effective
                double fg2 = v1 * v1 + v2 * v2;
                double r = Math.Sqrt(fg2);
                double rr = 1d / r;
                gp.c = Math.Abs(v1) * rr;
                gp.s = v2 * rr;

                if (v1 < 0) gp.s = -gp.s;
            }

            return gp;
        }

        /// <summary>
        /// Struct containing the cosine and sine elements of the Givens rotation.
        /// </summary>
        struct Givens
        {
            public double c;
            public double s;
        }
    }
}
