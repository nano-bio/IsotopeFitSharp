using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

using MathNet.Numerics.LinearAlgebra;

namespace IsotopeFit
{
    /// <summary>
    /// Static class containing definitions of all IsotopeFit data storage classes.
    /// </summary>
    public static class IFData
    {
        /// <summary>
        /// Utility function that converts a 2D array into a MathNet vector.
        /// </summary>
        /// <param name="arr">2D array that was read from an IFD file.</param>
        /// <returns>MathNet vector.</returns>
        internal static double[] Arr2DTo1D(double[][] arr)
        {
            bool dim = arr.Length > arr[0].Length;

            int length = dim ? arr.Length : arr[0].Length;

            double[] tmp = new double[length];

            for (int i = 0; i < length; i++)
            {
                if (dim)
                {
                    tmp[i] = arr[i][0];
                }
                else
                {
                    tmp[i] = arr[0][i];
                }
            }

            return tmp;
        }

        /// <summary>
        /// Converts a 2D array into a MathNet matrix.
        /// </summary>
        /// <param name="arr">2D array that was read from an IFD file.</param>
        /// <returns>MathNet matrix.</returns>
        internal static Matrix<double> Arr2DToMatrix(double[][] arr)
        {
            return Matrix<double>.Build.DenseOfRowArrays(arr);
        }

        /// <summary>
        /// Storage class for mass spectrum data.
        /// </summary>
        public class Spectrum
        {
            #region Fields

            private double[] rawMassAxis;
            private double[] rawSignalAxis;

            #endregion

            #region Properties

            public int RawLength { get; private set; }

            /// <summary>
            /// Raw experimental mass axis.
            /// </summary>
            public double[] RawMassAxis
            {
                get { return rawMassAxis; }
                set
                {
                    RawLength = value.Length;
                    rawMassAxis = value;
                }
            }

            /// <summary>
            /// Raw experimental signal axis.
            /// </summary>
            public double[] RawSignalAxis
            {
                get { return rawSignalAxis; }
                set
                {
                    RawLength = value.Length;
                    rawSignalAxis = value;
                }
            }

            /// <summary>
            /// Mass axis corrected for mass offset.
            /// </summary>
            public double[] MassAxis { get; internal set; }

            /// <summary>
            /// Signal axis with baseline subtracted.
            /// </summary>
            public double[] SignalAxis { get; internal set; }

            /// <summary>
            /// Crop start index.
            /// </summary>
            public int CropStartIndex { get; set; }

            /// <summary>
            /// Crop end index.
            /// </summary>
            public int CropEndIndex { get; set; }

            /// <summary>
            /// Crop start mass.
            /// </summary>
            public double CropStartMass { get; set; }

            /// <summary>
            /// Crop end mass.
            /// </summary>
            public double CropEndMass { get; set; }

            public double[] RawMassAxisCrop { get; set; }
            public double[] SignalAxisCrop { get; set; }
            public int CroppedLength { get; internal set; }
            internal bool Cropped { get; set; }

            #endregion

            #region Constructors

            /// <summary>
            /// Creates an empty Spectrum object, to be filled with data later.
            /// </summary>
            public Spectrum() { }

            /// <summary>
            /// Creates new storage for mass spectrum data from specified x and y axis.
            /// </summary>
            /// <param name="massAxis">Mass axis of the spectrum.</param>
            /// <param name="signalAxis">Signal axis of the spectrum.</param>
            public Spectrum(double[] massAxis, double[] signalAxis)
            {
                if (massAxis.Length != signalAxis.Length) throw new WorkspaceException("Supplied arrays have different lengths.");

                RawLength = massAxis.Length;
                RawMassAxis = massAxis;
                RawSignalAxis = signalAxis;
            }

            /// <summary>
            /// Creates new storage for raw data and populates it with supplied data.
            /// </summary>
            /// <remarks>The <paramref name="data"/> is required to have a Nx2 dimension.</remarks>
            /// <param name="data">Data to be stored in the new instance.</param>
            internal Spectrum(double[][] data)
            {
                if (data == null || data[0].Length != 2) throw new WorkspaceException("The supplied data are invalid.");

                RawLength = data.Length;
                RawMassAxis = new double[RawLength];
                RawSignalAxis = new double[RawLength];

                for (int i = 0; i < RawLength; i++)
                {
                    RawMassAxis[i] = data[i][0];
                    RawSignalAxis[i] = data[i][1];
                }
            }

            #endregion
        }

        /// <summary>
        /// Class for storing data about a single cluster/molecule/fragment.
        /// </summary>
        public class Cluster
        {
            /// <summary>
            /// Stores the isotopical data for the cluster.
            /// </summary>
            public IsotopeData PeakData { get; set; }

            /// <summary>
            /// Human readable name of the cluster. NOT the ID!
            /// </summary>
            public string Name { get; set; }

            //public double MinMass { get; set; }
            //public double MaxMass { get; set; }

            /// <summary>
            /// Centre of mass of the cluster, used for mass offset correction.
            /// </summary>
            /// <remarks>
            /// This value needs to be set only for those clusters, that were chosen as calibration clusters by the <see cref="Calibration.NameList"/> ID list.
            /// </remarks> 
            public double CentreOfMass { get; set; }

            //public int MinIndex { get; set; }
            //public int MaxIndex { get; set; }


            /// <summary>
            /// Fitted abundance of the cluster.
            /// </summary>
            public double Abundance { get; set; }

            /// <summary>
            /// Error of the fitted abundance value.
            /// </summary>
            public double AbundanceError { get; set; }

            /// <summary>
            /// Mass offset of the cluster, used for mass offset correction.
            /// </summary>
            /// <remarks>
            /// This value needs to be set only for those clusters, that were chosen as calibration clusters by the <see cref="Calibration.NameList"/> ID list.
            /// </remarks> 
            public double MassOffset { get; set; }

            /// <summary>
            /// Resolution of the particular cluster line, used for resolution calibration.
            /// </summary>
            /// <remarks>
            /// This value needs to be set only for those clusters, that were chosen as calibration clusters by the <see cref="Calibration.NameList"/> ID list.
            /// </remarks>
            public double Resolution { get; set; }

            //[Obsolete]
            //public int RootIndex { get; set; }

            /// <summary>
            /// Creates an empty cluster object.
            /// </summary>
            public Cluster()
            {
                PeakData = new IsotopeData();
            }

            /// <summary>
            /// Class to store data about the isotopical pattern of a molecule/fragment.
            /// </summary>
            public class IsotopeData
            {
                public double[] Mass { get; set; }
                public double[] Abundance { get; set; }

                public IsotopeData() { }

                public IsotopeData(double[][] data)
                {
                    Mass = new double[data.Length];
                    Abundance = new double[data.Length];

                    for (int i = 0; i < data.Length; i++)
                    {
                        Mass[i] = data[i][0];
                        Abundance[i] = data[i][1];
                    }
                }

                /// <summary>
                /// Method that adds a single isotopical entry to the <see cref="IsotopeData"/> object.
                /// </summary>
                /// <remarks>
                /// The input array is to contain 2-element array, where the first element is the mass and the second is the abundance.
                /// </remarks>
                /// <param name="input">Array of inputs.</param>
                public void Add(double[] input)
                {
                    // if the Mass and Abundance arrays create them and add the elements
                    if (Mass == null && Abundance == null)
                    {
                        Mass = new double[1] { input[0] };
                        Abundance = new double[1] { input[1] };
                    }
                    else    // if they exits, make a new longer arrays and append the new elements at their ends
                    {
                        double[] newMass = new double[Mass.Length + 1];
                        double[] newAbundance = new double[Mass.Length + 1];

                        Array.Copy(Mass, newMass, Mass.Length);
                        Array.Copy(Abundance, newAbundance, Abundance.Length);

                        newMass[Mass.Length] = input[0];
                        newAbundance[Abundance.Length] = input[0];

                        Mass = newMass;
                        Abundance = newAbundance;
                    }
                }
            }
        }

        public class Calibration
        {
            public double[] COMList { get; set; }
            public double[] MassOffsetList { get; set; }
            public double[] ResolutionList { get; set; }
            public string MassOffsetMethod { get; set; }
            public string ResolutionMethod { get; set; }
            public int MassOffsetParam { get; set; }
            public int ResolutionParam { get; set; }

            /// <summary>
            /// List of IDs of clusters picked as calibration points for mass offset and resolution.
            /// </summary>
            public List<string> NameList { get; set; }

            /// <summary>
            /// Object containing the information about the line shape.
            /// </summary>
            public LineShape Shape { get; set; }

            /// <summary>
            /// Object containing the results of mass offset interpolation.
            /// </summary>
            public Interpolation MassOffsetInterp { get; internal set; }

            /// <summary>
            /// Object containing the results of resolution fit.
            /// </summary>
            public Interpolation ResolutionInterp { get; internal set; }

            internal Calibration()
            {
                NameList = new List<string>();
                Shape = new LineShape();
            }

            /// <summary>
            /// Class containing the information about line shape.
            /// </summary>
            public class LineShape
            {
                public LineShape() { }

                public LineShape(double[] breaks, double[][] coeffs)
                {
                    Breaks = breaks;
                    Coeffs = Matrix<double>.Build.DenseOfRowArrays(coeffs);
                }

                //public string Form { get; set; }

                /// <summary>
                /// Breaks of the partial polynomial describing the line shape
                /// </summary>
                public double[] Breaks { get; set; }

                /// <summary>
                /// Coefficients of the partial polynomial describing the line shape
                /// </summary>
                public Matrix<double> Coeffs { get; set; }
                //public int Pieces { get; set; }
                //public int Order { get; set; }
                //public int Dim { get; set; }
            }
        }

        /// <summary>
        /// Class containing data for the baseline correction.
        /// </summary>
        public class BaselineCorr
        {
            internal BaselineCorr() { }

            //public double StartMass { get; set; }
            //public double EndMass { get; set; }
            public int NDiv { get; set; }
            public double Percent { get; set; }

            /// <summary>
            /// X-values of the baseline correction points.
            public double[] XAxis { get; set; }

            /// <summary>
            /// Y-values of the baseline correction points.
            /// </summary>
            public double[] YAxis { get; set; }

            /// <summary>
            /// Object containing results of the baseline interpolation.
            /// </summary>
            public Interpolation BaselineInterpolation { get; internal set; }            
        }

        /// <summary>
        /// Class containing additional information about the currently loaded file.
        /// </summary>
        public class FileInfo
        {
            internal FileInfo() { }

            /// <summary>
            /// Name of the H5 file that contains the measurement data.
            /// </summary>
            public string OriginalFilename { get; set; }

            /// <summary>
            /// Full path to the H5 file that contains the measurement data.
            /// </summary>
            public string H5CompleteName { get; set; }            
        }
    }
}
