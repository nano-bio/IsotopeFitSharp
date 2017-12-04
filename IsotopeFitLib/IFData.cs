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
        [Obsolete]
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
        [Obsolete]
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
            public int CroppedLength { get; set; }
            //TODO: crop start and crop end?

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
            /// Crop start index. Not implemented at the moment.
            /// </summary>
            public int CropStartIndex { get; set; }

            /// <summary>
            /// Crop end index. Not implemented at the moment.
            /// </summary>
            public int CropEndIndex { get; set; }

            ///// <summary>
            ///// Contains the last calculated baseline that was subtracted from <see cref="Spectrum.RawSignalAxis"/> to obtain <see cref="Spectrum.SignalAxis"/>.
            ///// </summary>
            //public double[] Baseline { get; internal set; }

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
            /// Name of the cluster. NOT the ID!
            /// </summary>
            public string Name { get; set; }

            public double MinMass { get; set; }
            public double MaxMass { get; set; }

            /// <summary>
            /// Centre of mass of the cluster, used for mass offset correction.
            /// </summary>
            /// <remarks>
            /// This value needs to be set only for those clusters, that were chosen as calibration clusters by the <see cref="Calibration.NameList"/> ID list.
            /// </remarks> 
            public double CentreOfMass { get; set; }

            public int MinIndex { get; set; }
            public int MaxIndex { get; set; }

            public double Area { get; set; }
            public double AreaError { get; set; }

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

            [Obsolete]
            public int RootIndex { get; set; }

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

                internal IsotopeData() { }

                internal IsotopeData(double[][] data)
                {
                    Mass = new double[data.Length];
                    Abundance = new double[data.Length];

                    for (int i = 0; i < data.Length; i++)
                    {
                        Mass[i] = data[i][0];
                        Abundance[i] = data[i][1];
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
            public List<string> NameList { get; set; }
            public LineShape Shape { get; set; }

            public Interpolation MassOffsetInterp { get; internal set; }
            public Interpolation ResolutionInterp { get; internal set; }

            internal Calibration()
            {
                NameList = new List<string>();
                Shape = new LineShape();
            }

            public class LineShape
            {
                public string Form { get; set; }
                public double[] Breaks { get; set; }
                public Matrix<double> Coefs { get; set; }
                public int Pieces { get; set; }
                public int Order { get; set; }
                public int Dim { get; set; }
            }
        }

        public class BaselineCorr
        {
            public double StartMass { get; set; }
            public double EndMass { get; set; }
            public int NDiv { get; set; }
            public double Percent { get; set; }
            public double[] XAxis { get; set; }
            public double[] YAxis { get; set; }

            public Interpolation BaselineInterpolation { get; internal set; }

            internal BaselineCorr()
            {

            }
        }
    }
}
