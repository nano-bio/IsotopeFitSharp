using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Storage;

//using CSparse;
using CSparse.Double;
using CSparse.Storage;

namespace IsotopeFit
{
    public static partial class Algorithm
    {
        /// <summary>
        /// Calculates QR factorization of a least squares system using Householder transformation.
        /// </summary>
        /// <remarks>
        /// Uses QR factorization implementation from CSparse library.
        /// </remarks>
        /// <param name="lss">Least squares system object to be QR factorized.</param>
        /// <returns>QR factorized least squares system object.</returns>        
        public static LeastSquaresSystem LeaSqrSparseQRHouseholder(LeastSquaresSystem lss)
        {
            //TODO: cleanup 
            //SparseMatrix M = (SparseMatrix)SparseMatrix.Build.SparseOfMatrix(lss.DesignMatrix);
            CSparse.Double.SparseMatrix M;

            if (lss.ObservationVector != null)
            {
                SparseVector v = (SparseVector)SparseVector.Build.SparseOfVector(lss.ObservationVector);    //TODO: this is also a little wasteful
                SparseVectorStorage<double> svs = v.Storage as SparseVectorStorage<double>;

                //TODO: and we have to put the vector as a last column of the matrix to be factorized

                int nonZeroCount = lss.DesignMatrix.NonZerosCount + svs.ValueCount;

                double[] values = new double[nonZeroCount];
                int[] rowIndices = new int[nonZeroCount];
                int[] colPointers = new int[lss.DesignMatrix.ColumnCount + 2];    // +1 because n-th colPointer tells about the number of elements in (n-1)st column

                Array.Copy(lss.DesignMatrix.Values, values, lss.DesignMatrix.Values.Length);
                Array.Copy(svs.Values, 0, values, lss.DesignMatrix.Values.Length, svs.Values.Length);

                Array.Copy(lss.DesignMatrix.RowIndices, rowIndices, lss.DesignMatrix.RowIndices.Length);
                Array.Copy(svs.Indices, 0, rowIndices, lss.DesignMatrix.RowIndices.Length, svs.Indices.Length);

                Array.Copy(lss.DesignMatrix.ColumnPointers, colPointers, lss.DesignMatrix.ColumnPointers.Length);
                colPointers[colPointers.Length - 1] = colPointers[colPointers.Length - 2] + svs.ValueCount;

                M = new CSparse.Double.SparseMatrix(lss.DesignMatrix.RowCount, lss.DesignMatrix.ColumnCount + 1)
                {
                    Values = values,
                    RowIndices = rowIndices,
                    ColumnPointers = colPointers
                };
            }
            else // otherwise it means the column has been added (it is the caller responsibility that the data make sense)
            {
                M = lss.DesignMatrix;
            }


            // not necessary, because we already have the CSparse CCS format
            //M = (SparseMatrix)M.Transpose();
            //SparseCompressedRowMatrixStorage<double> spStor = (SparseCompressedRowMatrixStorage<double>)M.Storage;

            //var C = new CSparse.Double.SparseMatrix(M.ColumnCount, M.RowCount)
            //{
            //    ColumnPointers = spStor.RowPointers,
            //    RowIndices = spStor.ColumnIndices,
            //    Values = spStor.Values
            //};

            var QR = CSparse.Double.Factorization.SparseQR.Create(M, CSparse.ColumnOrdering.MinimumDegreeAtA);
            var R = (CSparse.Double.SparseMatrix)QR.R;

            R.DropZeros();  //TODO: this might need to be set to machine epsilon

            // now we need to cut the solution
            // access by matrix indices does not work with the CSparse format, we have to copy the internal arrays
            CSparse.Double.SparseMatrix Mr = new CSparse.Double.SparseMatrix(R.ColumnCount - 1, R.ColumnCount - 1)
            {
                Values = new double[R.ColumnPointers[R.ColumnPointers.Length - 2]],
                RowIndices = new int[R.ColumnPointers[R.ColumnPointers.Length - 2]],
                ColumnPointers = new int[R.ColumnPointers.Length - 1]
            };

            Array.Copy(R.Values, Mr.Values, Mr.Values.Length);
            Array.Copy(R.RowIndices, Mr.RowIndices, Mr.RowIndices.Length);
            Array.Copy(R.ColumnPointers, Mr.ColumnPointers, Mr.ColumnPointers.Length);

            //for (int i = 0; i < Mr.RowCount; i++)
            //{
            //    for (int j = 0; j < Mr.ColumnCount; j++)
            //    {
            //        //Mr.At(i, j, R.At(i, j));
            //        Mr.ke
            //    }
            //}

            double[] kurva = new double[Mr.ColumnCount];
            Array.Copy(R.Column(R.ColumnCount - 1), 0, kurva, 0, Mr.ColumnCount);
            Vector<double> vr = Vector<double>.Build.DenseOfArray(kurva);

            LeastSquaresSystem factorizedLSS = new LeastSquaresSystem(Mr, vr)
            {
                ColumnOrdering = new int[QR.columnOrdering.Length]
            };

            Array.Copy(QR.columnOrdering, factorizedLSS.ColumnOrdering, QR.columnOrdering.Length);

            //QR = null;

            return factorizedLSS;
        }

        ///// <summary>
        ///// Calculates QR factorization of a least squares system using Givens rotations.
        ///// </summary>
        ///// <param name="lss">Least squares system object to be QR factorized.</param>
        ///// <returns>QR factorized least squares system object.</returns>
        //public static LeastSquaresSystem LeaSqrSparseQRGivens(LeastSquaresSystem lss)
        //{
        //    Matrix<double> M = Matrix<double>.Build.SparseOfMatrix(lss.DesignMatrix);
        //    Vector<double> v = Vector<double>.Build.DenseOfVector(lss.ObservationVector);

        //    Givens gp;
        //    double x, y;

        //    for (int i = 0; i < M.ColumnCount; i++)
        //    {
        //        for (int j = M.RowCount - 1; j > i; j--)
        //        {
        //            if (M.At(j, i) == 0) continue;

        //            gp = GivensParams(M.At(i, i), M.At(j, i));

        //            for (int k = i; k < M.ColumnCount; k++)
        //            {
        //                x = M.At(i, k);
        //                y = M.At(j, k);

        //                M.At(i, k, gp.c * x + gp.s * y);
        //                M.At(j, k, gp.c * y - gp.s * x);

        //                //TODO: add zero-coercing
        //            }

        //            x = v.At(i);
        //            y = v.At(j);

        //            //TODO: this must be tested, if it is done correctly
        //            v.At(i, gp.c * x + gp.s * y);
        //            v.At(j, gp.c * y - gp.s * x);
        //        }
        //    }

        //    //TODO: this cutting, especially for the vector should maybe be checked for non-zero values before the cutting line
        //    return new LeastSquaresSystem(M.SubMatrix(0, M.ColumnCount, 0, M.ColumnCount), v.SubVector(0, M.ColumnCount));
        //}

        ///// <summary>
        ///// Calculates QR factorization of the sparse matrix M using Householder transformation.
        ///// </summary>
        ///// <remarks>
        ///// Method calculates the upper triangular matrix R without explicitly forming the Q matrix.
        ///// It is meant primarily for sparse matrices.
        ///// </remarks>
        ///// <param name="M">Matrix to be QR factorized.</param>
        //public static Matrix<double> SparseQRHouseholder(Matrix<double> M)
        //{
        //    M = (SparseMatrix)M.Transpose();
        //    SparseCompressedRowMatrixStorage<double> spStor = (SparseCompressedRowMatrixStorage<double>)M.Storage;

        //    var C = new CSparse.Double.SparseMatrix(M.ColumnCount, M.RowCount)
        //    {
        //        ColumnPointers = spStor.RowPointers,
        //        RowIndices = spStor.ColumnIndices,
        //        Values = spStor.Values
        //    };

        //    var R = CSparse.Double.Factorization.SparseQR.Create(C, CSparse.ColumnOrdering.MinimumDegreeAtA).R;
        //    R.DropZeros();  //TODO: this might need to be set to machine epsilon

        //    Matrix<double> Mr = Matrix<double>.Build.Dense(R.ColumnCount, R.ColumnCount);

        //    for (int i = 0; i < Mr.RowCount; i++)
        //    {
        //        for (int j = 0; j < Mr.ColumnCount; j++)
        //        {
        //            Mr.At(i, j, R.At(i, j));
        //        }
        //    }

        //    return Mr;
        //}

        /// <summary>
        /// Calculates QR factorization of the sparse matrix M using Givens rotations.
        /// </summary>
        /// <remarks>
        /// Method calculates the upper triangular matrix R without explicitly forming the Q matrix.
        /// It is meant primarily for sparse matrices.
        /// </remarks>
        /// <param name="M">Matrix to be QR factorized.</param>
        public static Matrix<double> SparseQRGivens(Matrix<double> M)
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
                //TODO: test if putting these operations together is more effective
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
