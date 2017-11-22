using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

using MathNet.Numerics.LinearAlgebra;

using IsotopeFit.Numerics;

namespace IsotopeFit
{
    public partial class Workspace
    {
        internal class DesignMtrx
        {
            #region Fields

            private Vector<double>[] designMatrixVectors;
            private double[] massAxis;

            bool[] fitMask;

            PolyInterpolation resolutionFit;

            #endregion

            #region Constructors

            /// <summary>
            /// Create new design matrix and initialize data sources for further calculations.
            /// </summary>
            internal DesignMtrx(IFData.Spectrum spectrum, List<IFData.Molecule> molecules, IFData.Calibration calibration)
            {
                massAxis = spectrum.MassAxis.ToArray();
                Molecules = molecules;
                Calibration = calibration;

                Rows = spectrum.Length;
                Cols = molecules.Count;
            }

            #endregion

            #region Properties

            //private Vector<double> MassAxis { get; set; }
            private List<IFData.Molecule> Molecules { get; set; }
            private IFData.Calibration Calibration { get; set; }            

            internal Matrix<double> Storage { get; private set; }   //TODO: maybe a field would suffice and change it directly to a sparse matrix
            internal int Rows { get; private set; }
            internal int Cols { get; private set; }
            internal Matrix<double> R { get; private set; }

            internal double SearchRange { get; set; }   //TODO: at the moment, those two values are not being set anywhere
            internal double FwhmRange { get; set; }

            #endregion

            //private delegate Vector<double> BuildColumnDelegate(int index);
            //private Action<int> buildColumnDelegate;

            #region Methods

            internal void Build()
            {                
                //TODO: make them sparse
                designMatrixVectors = new Vector<double>[Cols];
                fitMask = new bool[Rows];

                //TODO: this is temporary
                SearchRange = 1;
                FwhmRange = 0.5;

                // initialize the resolution interpolation object from IFD data. immediately calculates the interpolation internally
                resolutionFit = new PolyInterpolation(Calibration.COMList.ToArray(), Calibration.ResolutionList.ToArray(), (int)Calibration.ResolutionParam); //TODO: maybe change the ulongs to ints

                // loop through all molecules
                Parallel.For(0, Cols, BuildInit, BuildWork, BuildFinal);

                //build the sparse design matrix from the column vector array
                Storage = Matrix<double>.Build.SparseOfColumnVectors(designMatrixVectors);   //TODO: this is horribly slow. we should build the matrix manually.

            }


            /// <summary>
            /// Initializes a workspace for each thread of the parallel for loop that builds the design matrix.
            /// Executes once, before the first iteration runs on a newly created thread.
            /// </summary>
            /// <returns>New instance of thread workspace.</returns>
            private BuildState BuildInit()
            {                
                return new BuildState(Rows, resolutionFit.Coefs);
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
                int isotopePeakCount = Molecules[moleculeIndex].PeakData.Mass.Count;

                // loop through all isotope peaks of the current molecule
                for (int i = 0; i < isotopePeakCount; i++)
                {
                    int peakLowerLimitIndex, peakUpperLimitIndex;
                    int fitmaskLowerLimitIndex, fitmaskUpperLimitIndex;
                    double mass, abundance, resolution, fwhm;

                    Vector<double> breaks;
                    Matrix<double> coefs;

                    Vector<double> currentColumn;

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
                    currentColumn = Vector<double>.Build.Sparse(Rows);  //TODO: maybe we can make this building also manually, to be more effective
                    PPInterpolation peakshape = new PPInterpolation(breaks.ToArray(), coefs.ToRowArrays());
                    
                    for (int j = peakLowerLimitIndex; j <= peakUpperLimitIndex; j++)
                    {
                        // check if the point falls within the fitmask
                        if (bs.fitMask[j])
                        {
                            // evaluate the peakshape partial polynomial at the j index of mass axis and add the point to the design matrix column
                            currentColumn[j] = abundance * peakshape.Evaluate(massAxis[j]);
                        }
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
                Vector<double> breaks = Vector<double>.Build.Dense(sh.Breaks.Count, 0); 

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

                for (int row = 0; row < coefs.RowCount; row++)
                {
                    for (int col = 0; col < coefs.ColumnCount; col++)
                    {
                        // value 4 is hardcoded, because the line shape is always represented by piecewise cubic polynomials
                        switch (col % 4)
                        {
                            case 0:
                                coefs.At(row, col, sh.Coefs.At(row, col) * abundance / fwhm);
                                break;
                            case 1:
                                coefs.At(row, col, sh.Coefs.At(row, col) * abundance / Math.Pow(fwhm, 2));
                                break;
                            case 2:
                                coefs.At(row, col, sh.Coefs.At(row, col) * abundance / Math.Pow(fwhm, 3));
                                break;
                            case 3:
                                coefs.At(row, col, sh.Coefs.At(row, col) * abundance / Math.Pow(fwhm, 4));
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
                internal PolyInterpolation resolutionFit;

                internal BuildState(int size, double[] resFitCoefs)
                {
                    fitMask = new bool[size];
                    resolutionFit = new PolyInterpolation(resFitCoefs);
                }
            }

            #endregion
        }
    }
}
