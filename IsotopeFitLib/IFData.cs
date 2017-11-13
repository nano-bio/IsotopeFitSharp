using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet.Numerics.LinearAlgebra;

namespace IsotopeFit
{
    public static class IFData
    {
        //TODO: those vectors will maybe have to be changed to C arrays for better python compatibility

        /// <summary>
        /// Utility function that converts a 2D array into a MathNet vector.
        /// </summary>
        /// <param name="arr">2D array that was read from an IFD file.</param>
        /// <returns>MathNet vector.</returns>
        [Obsolete]
        internal static Vector<double> Arr2DToVect(double[][] arr)
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

            return Vector<double>.Build.DenseOfArray(tmp);
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
            public int Length { get; set; }

            public Vector<double> MassAxis { get; set; }
            public Vector<double> SignalAxis { get; set; }

            internal Spectrum()
            {

            }

            /// <summary>
            /// Creates new storage for mass spectrum data from specified x and y axis.
            /// </summary>
            /// <param name="massAxis">Mass axis of the spectrum.</param>
            /// <param name="signalAxis">Signal axis of the spectrum.</param>
            internal Spectrum(Vector<double> massAxis, Vector<double> signalAxis)
            {
                Length = massAxis.Count;
                MassAxis = massAxis;
                SignalAxis = signalAxis;
            }

            /// <summary>
            /// Creates new storage for raw data and populates it with supplied data.
            /// </summary>
            /// <param name="data">Data to be stored in the new instance.</param>
            internal Spectrum(double[][] data)
            {
                //TODO: include the check to determine is the array is MxN or NxM

                Length = data.Length;
                MassAxis = Vector<double>.Build.Dense(Length, 0);
                SignalAxis = Vector<double>.Build.Dense(Length, 0);

                for (int i = 0; i < Length; i++)
                {
                    MassAxis[i] = data[i][0];
                    SignalAxis[i] = data[i][1];
                }
            }
        }

        /// <summary>
        /// Class for storing data about a single molecule/fragment.
        /// </summary>
        public class Molecule
        {
            public IsotopeData PeakData { get; set; }
            public string Name { get; set; }
            public double MinMass { get; set; }
            public double MaxMass { get; set; }
            public double CentreOfMass { get; set; }
            public ulong MinIndex { get; set; }
            public ulong MaxIndex { get; set; }
            public double Area { get; set; }
            public double AreaError { get; set; }
            public ulong RootIndex { get; set; }

            internal Molecule()
            {
                //PeakData = new IsotopeData();
            }

            /// <summary>
            /// Class to store data about the isotopical pattern of a molecule/fragment.
            /// </summary>
            public class IsotopeData
            {
                public Vector<double> Mass { get; set; }
                public Vector<double> Abundance { get; set; }

                internal IsotopeData(double[][] data)
                {
                    Mass = Vector<double>.Build.Dense(data.Length, 0);
                    Abundance = Vector<double>.Build.Dense(data.Length, 0);

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
            public Vector<double> COMList { get; set; }
            public Vector<double> MassOffsetList { get; set; }
            public Vector<double> ResolutionList { get; set; }
            public string MassOffsetMethod { get; set; }
            public string ResolutionMethod { get; set; }
            public ulong MassOffsetParam { get; set; }
            public ulong ResolutionParam { get; set; }
            public List<string> Namelist { get; set; }
            public LineShape Shape { get; set; }

            internal Calibration()
            {
                Namelist = new List<string>();
            }

            public class LineShape
            {
                public string Form { get; set; }
                public Vector<double> Breaks { get; set; }
                public Matrix<double> Coefs { get; set; }
                public ulong Pieces { get; set; }
                public ulong Order { get; set; }
                public ulong Dim { get; set; }
            }
        }

        public class BaselineCorr
        {
            public double StartMass { get; set; }
            public double EndMass { get; set; }
            public ulong NDiv { get; set; }
            public double Percent { get; set; }
            public Vector<double> XAxis { get; set; }
            public Vector<double> YAxis { get; set; }

            internal BaselineCorr()
            {

            }
        }
    }
}
