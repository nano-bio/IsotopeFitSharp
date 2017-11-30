using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Storage;

using CSparse.Double;

using IsotopeFit.Numerics;

namespace IsotopeFit
{
    public partial class Workspace
    {
        public class DesignMtrx
        {
            #region Fields

            private MathNet.Numerics.LinearAlgebra.Double.SparseVector[] designMatrixVectors;
            private MathNet.Numerics.LinearAlgebra.Double.SparseVector observationVector;
            private double[] massAxis;

            bool[] fitMask;

            Interpolation resolutionFit;

            #endregion

            #region Constructors

            /// <summary>
            /// Create new design matrix and initialize data sources for further calculations.
            /// </summary>
            internal DesignMtrx(IFData.Spectrum spectrum, List<IFData.Molecule> molecules, IFData.Calibration calibration, Interpolation resolutionInterp)
            {
                massAxis = spectrum.MassAxis.ToArray();
                observationVector = (MathNet.Numerics.LinearAlgebra.Double.SparseVector)MathNet.Numerics.LinearAlgebra.Double.SparseVector.Build.SparseOfArray(spectrum.SignalAxis);
                Molecules = molecules;
                Calibration = calibration;
                resolutionFit = resolutionInterp;

                Rows = spectrum.RawLength;
                Cols = molecules.Count;
            }

            #endregion

            #region Properties

            //private Vector<double> MassAxis { get; set; }
            private List<IFData.Molecule> Molecules { get; set; }
            private IFData.Calibration Calibration { get; set; }            

            public SparseMatrix Storage { get; private set; }   //TODO: maybe a field would suffice and change it directly to a sparse matrix
            internal int Rows { get; private set; }
            internal int Cols { get; private set; }
            internal Matrix<double> R { get; private set; }

            internal double SearchRange { get; set; }   //TODO: at the moment, those two values are not being set anywhere
            internal double FwhmRange { get; set; }

            #endregion

            #region Methods

            internal void Build()
            {
                //TODO: make them sparse
                designMatrixVectors = new MathNet.Numerics.LinearAlgebra.Double.SparseVector[Cols + 1]; // plus 1, because that will be the observations vector. this makes later calculations faster.
                fitMask = new bool[Rows];

                //TODO: this is temporary
                SearchRange = 1;
                FwhmRange = 0.5;

                // loop through all molecules
                Parallel.For(0, Cols, BuildInit, BuildWork, BuildFinal);


                // add the observation vector to the last column
                int fitMaskNonZero = 0;
                for (int i = 0; i < fitMask.Length; i++)
                {
                    if (fitMask[i]) fitMaskNonZero++;
                }
                
                double[] values2 = new double[fitMaskNonZero];
                int[] indices2 = new int[fitMaskNonZero];

                int j = 0;
                for (int i = 0; i < fitMask.Length; i++)
                {
                    if (fitMask[i])
                    {
                        values2[j] = (observationVector.Storage as SparseVectorStorage<double>).Values[i];
                        indices2[j] = (observationVector.Storage as SparseVectorStorage<double>).Indices[i];
                        j++;
                    }
                }

                SparseVectorStorage<double> obsVecStor = SparseVectorStorage<double>.OfValue(fitMaskNonZero, 1);
                obsVecStor.Values = values2;
                obsVecStor.Indices = indices2;

                designMatrixVectors[Cols] = (MathNet.Numerics.LinearAlgebra.Double.SparseVector)MathNet.Numerics.LinearAlgebra.Double.SparseVector.Build.Sparse(obsVecStor);
                //designMatrixVectors[Cols] = (MathNet.Numerics.LinearAlgebra.Double.SparseVector)MathNet.Numerics.LinearAlgebra.Double.SparseVector.Build.Sparse(10);

                //build the sparse design matrix from the column vector array                
                //TODO: since we are building the matrix from columns, it should be easy to build a compressed-column-storage right away. That is the storage type needed for the QR factorization.
                //Storage = Matrix<double>.Build.SparseOfColumnVectors(designMatrixVectors);
                
                // fitmask application
                for (int i = 0; i < designMatrixVectors.Length - 1; i++)
                {                    
                    SparseVectorStorage<double> st = designMatrixVectors[i].Storage as SparseVectorStorage<double>;

                    Array.Resize(ref st.Indices, st.ValueCount);    // TODO: there are useless trailing zeroes, we need to cut them out

                    int minIndex = st.Indices[0];
                    int length = st.Indices.Length;

                    List<int> nonMaskedIndices = new List<int>(length);
                    List<double> nonMaskedValues = new List<double>(length);

                    // extract the masked indices and values (strictly speaking, the complement to masked indices)
                    for (int k = 0; k < length; k++)
                    {
                        if (fitMask[st.Indices[k]])
                        {
                            nonMaskedIndices.Add(st.Indices[k]);
                            nonMaskedValues.Add(st.Values[k]);
                        }
                    }

                    // put the masked indices and values in the storage object, replacing the non-masked data
                    (designMatrixVectors[i].Storage as SparseVectorStorage<double>).Indices = nonMaskedIndices.ToArray();
                    (designMatrixVectors[i].Storage as SparseVectorStorage<double>).Values = nonMaskedValues.ToArray();
                    (designMatrixVectors[i].Storage as SparseVectorStorage<double>).ValueCount = nonMaskedIndices.Count;
                }

                // create the design matrix
                int nonZeroCount = 0;

                for (int i = 0; i < designMatrixVectors.Length; i++)
                {
                    nonZeroCount += designMatrixVectors[i].NonZerosCount;
                }

                double[] values = new double[nonZeroCount];
                int[] rowIndices = new int[nonZeroCount];
                int[] colPointers = new int[designMatrixVectors.Length + 1];    // +1 because n-th colPointer tells about the number of elements in (n-1)st column

                for (int i = 0; i < designMatrixVectors.Length; i++)
                {
                    Array.Copy((designMatrixVectors[i].Storage as SparseVectorStorage<double>).Values, 0, values, colPointers[i], designMatrixVectors[i].NonZerosCount);
                    Array.Copy((designMatrixVectors[i].Storage as SparseVectorStorage<double>).Indices, 0, rowIndices, colPointers[i], designMatrixVectors[i].NonZerosCount);
                    colPointers[i + 1] = colPointers[i] + designMatrixVectors[i].NonZerosCount;
                }

                Storage = new SparseMatrix(Rows, Cols + 1)
                {
                    Values = values,
                    RowIndices = rowIndices,
                    ColumnPointers = colPointers
                };
            }


            /// <summary>
            /// Initializes a workspace for each thread of the parallel for loop that builds the design matrix.
            /// Executes once, before the first iteration runs on a newly created thread.
            /// </summary>
            /// <returns>New instance of thread workspace.</returns>
            private BuildState BuildInit()
            {                
                return new BuildState(Rows, resolutionFit);
            }

            /// <summary>
            /// Builds one column of the design matrix, which corresponds to a particular cluster.
            /// Executes once per iteration.
            /// </summary>
            /// <param name="moleculeIndex">Parallel for loop interation variable.</param>
            /// <param name="pls">Intra-thread messaging structure.</param>
            /// <param name="bs">Thread-local status object.</param>
            /// <returns>Modified thread status object.</returns>
            private BuildState BuildWork(int moleculeIndex, ParallelLoopState pls, BuildState bs)
            {
                //TODO: maybe we could generate only indices and values and return just that. it gets assebled to a matrix manually anyway. might have lower ram demands.

                int isotopePeakCount = Molecules[moleculeIndex].PeakData.Mass.Length;
                MathNet.Numerics.LinearAlgebra.Double.SparseVector currentColumn = (MathNet.Numerics.LinearAlgebra.Double.SparseVector)MathNet.Numerics.LinearAlgebra.Double.SparseVector.Build.Sparse(Rows);  //TODO: we can make this building also manually, to be more effective - and better as well

                // loop through all isotope peaks of the current molecule
                for (int i = 0; i < isotopePeakCount; i++)
                {
                    int peakLowerLimitIndex, peakUpperLimitIndex;
                    int fitmaskLowerLimitIndex, fitmaskUpperLimitIndex;
                    double mass, abundance, resolution, fwhm;

                    Vector<double> breaks;
                    Matrix<double> coefs;

                    mass = Molecules[moleculeIndex].PeakData.Mass[i];
                    abundance = Molecules[moleculeIndex].PeakData.Abundance[i];  // area of the line and abundance are proportional
                    resolution = bs.resolutionFit.Evaluate(mass);
                    fwhm = mass / resolution;

                    if (fwhm <= 0)
                    {
                        throw new ArithmeticException("fwhm was negative"); //TODO: reorder breaks and coefs? check matlab code
                    }

                    // transform shape breaks and coefficients
                    breaks = TransformLineShapeBreaks(Calibration.Shape, fwhm, mass);
                    coefs = TransformLineShapeCoefs(Calibration.Shape, fwhm, abundance);

                    // find peak and fitmask limits
                    peakLowerLimitIndex = FindLowerLimitIndex(massAxis, breaks.First());
                    peakUpperLimitIndex = FindUpperLimitIndex(massAxis, breaks.Last());
                    fitmaskLowerLimitIndex = FindLowerLimitIndex(massAxis, mass - SearchRange * FwhmRange * fwhm);
                    fitmaskUpperLimitIndex = FindUpperLimitIndex(massAxis, mass + SearchRange * FwhmRange * fwhm);

                    // adjust local fitmask
                    for (int j = fitmaskLowerLimitIndex; j <= fitmaskUpperLimitIndex; j++)
                    {
                        bs.fitMask[j] = true;
                    }

                    // build the design matrix column for the current molecule                    
                    PPInterpolation peakshape = new PPInterpolation(breaks.ToArray(), coefs.ToRowArrays());
                    
                    for (int j = peakLowerLimitIndex; j <= peakUpperLimitIndex; j++)
                    {
                        // TODO: WRONG, i have to calculate all of the points, because they might get unlocked by other fragments fitmasks. fitmask has to be applied only in the end
                        //if (bs.fitMask[j])
                        //{
                            // evaluate the peakshape partial polynomial at the j index of mass axis and add the point to the design matrix column
                            currentColumn[j] += abundance * peakshape.Evaluate(massAxis[j]);
                        //}
                    }

                    // add the built column to the column storage
                    designMatrixVectors[moleculeIndex] = currentColumn;
                }

                return bs;
            }

            /// <summary>
            /// Incorporates thread calculation results into the class-level storage.
            /// Executes once, after the last iteration of a particular thread.
            /// </summary>
            /// <param name="bs">Thread-local state object.</param>
            private void BuildFinal(BuildState bs)
            {
                lock (fitMask)
                {
                    for (int i = 0; i < Rows; i++)
                    {
                        if (bs.fitMask[i]) fitMask[i] = true;
                    }
                }
            }

            /// <summary>
            /// Calculates breaks of the piecewise polynomials that describe the line shape, adjusted for actual fwhm and mass.
            /// </summary>
            /// <param name="sh">Object containing information about the line shape.</param>
            /// <param name="fwhm">Value of full width at half maximum for the particular line, shape of which is being calculated.</param>
            /// <param name="mass">Mass of the fragment, to which the line being calculated corresponds.</param>
            /// <returns>Mathnet vector of recalculated piecewise polynomial breaks.</returns>
            private Vector<double> TransformLineShapeBreaks(IFData.Calibration.LineShape sh, double fwhm, double mass)
            {
                Vector<double> breaks = Vector<double>.Build.Dense(sh.Breaks.Length, 0); 

                for (int i = 0; i < breaks.Count; i++)
                {
                    breaks[i] = sh.Breaks[i] * fwhm + mass;
                }

                return breaks;
            }

            /// <summary>
            /// Calculates coefficients of the line shape piecewise polynomials, adjusted for actual fwhm and line area.
            /// </summary>
            /// <param name="sh">Object containing information about the line shape.</param>
            /// <param name="fwhm">Value of full width at half maximum for the particular line, shape of which is being calculated.</param>
            /// <param name="abundance">Isotopical abundance of the particular fragment, line for which is being calculated.</param>
            /// <returns>Mathnet matrix of recalculated piecewise polynomial coefficients.</returns>
            private Matrix<double> TransformLineShapeCoefs(IFData.Calibration.LineShape sh, double fwhm, double abundance)
            {
                Matrix<double> coefs = Matrix<double>.Build.Dense(sh.Coefs.RowCount, sh.Coefs.ColumnCount, 0);

                abundance = 1;  //TODO: remove

                for (int row = 0; row < coefs.RowCount; row++)
                {
                    for (int col = 0; col < coefs.ColumnCount; col++)
                    {
                        // value 4 is hardcoded, because the line shape is always represented by piecewise cubic polynomials
                        switch (col % 4)
                        {
                            case 0:
                                coefs.At(row, col, (sh.Coefs.At(row, col) * abundance / fwhm));  //fwhmDec
                                break;
                            case 1:
                                coefs.At(row, col, (sh.Coefs.At(row, col) * abundance / (fwhm * fwhm)));   //Math.Pow(fwhm, 2)
                                break;
                            case 2:
                                coefs.At(row, col, (sh.Coefs.At(row, col) * abundance / (fwhm * fwhm * fwhm)));   //Math.Pow(fwhm, 3)
                                break;
                            case 3:
                                coefs.At(row, col, (sh.Coefs.At(row, col) * abundance / (fwhm * fwhm * fwhm * fwhm)));   //Math.Pow(fwhm, 4)
                                break;
                        }
                    }
                }

                return coefs;    
            }

            /// <summary>
            /// Finds the index of the first value in a sorted array that is greater than the specified threshold.
            /// </summary>
            /// <param name="array">Sorted array of values to be searched.</param>
            /// <param name="threshold">Value to which the comparison will be made.</param>
            /// <returns>Index of first value that is greater than the specified threshold.</returns>
            private int FindLowerLimitIndex(double[] array, double threshold)
            {
                int index = Array.BinarySearch(array, threshold);

                if (index >= 0) return index;
                else return ~index;
            }

            /// <summary>
            /// Finds the index of the last value in a sorted array that is less than the specified threshold.
            /// </summary>
            /// <param name="array">Sorted array of values to be searched.</param>
            /// <param name="threshold">Value to which the comparison will be made.</param>
            /// <returns>Index of last value that is less than the specified threshold.</returns>
            private int FindUpperLimitIndex(double[] array, double threshold)
            {
                int index = Array.BinarySearch(array, threshold);

                if (index >= 0) return index;
                else return ~index - 1;  // -1 because we want the inner bound of the interval                
            }

            internal void CalculateQR()
            {
                throw new NotImplementedException();
            }

            #endregion

            #region Classes

            /// <summary>
            /// Class that contains iteration-specific and therefore thread-specific variables.
            /// </summary>
            private class BuildState
            {
                internal bool[] fitMask;
                internal Interpolation resolutionFit;
                internal int nonZeroCount;  //TODO: this might be faster solution, but lets leave it for now

                internal BuildState(int size, Interpolation resFit)
                {
                    fitMask = new bool[size];

                    switch (resFit.GetType().Name)
                    {
                        case "PolyInterpolation":
                            resolutionFit = resFit as PolyInterpolation;
                            break;
                        case "PPInterpolation":
                            resolutionFit = resFit as PPInterpolation;
                            break;
                    }
                    
                    nonZeroCount = 0;
                }
            }

            #endregion
        }
    }
}
