using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

using IsotopeFit.Numerics;

namespace IsotopeFit
{
   /*
    * TODO:
    * 1) zbastlit design matrix a ako posledny stlpec dat vektor pozorovani
    * 2) transponovat a napchat do csparse
    * 3) urobit QR
    * 4) vyrypat R a narezat ho na skutocne R a faktorizovany vektor pozorovani
    * 5) pokracovat v NNLS
    * 
    * Volitelne skusit ako si s ulohou poradi CSparse leasqr solver.
    */

    [ComVisible(true)]  //TODO: the ComVisible attribute can be set globally for the whole library. might be nicer.
    public partial class Workspace
    {
        #region Fields
        private DesignMtrx designMatrix;

        #endregion

        #region Constructors

        /// <summary>
        /// Create an empty IsotopeFit workspace.
        /// </summary>
        public Workspace()
        {
            //throw new NotImplementedException();
        }

        /// <summary>
        /// Create an IsotopeFit workspace and load an IFD file into it.
        /// </summary>
        /// <param name="IFDfile">Path to the IFD file to be loaded.</param>
        public Workspace(string path)
        {
            LoadIFDFile(path);
        }

        #endregion

        #region Destructors

        #endregion

        #region Properties

        public IFData.Spectrum SpectralData { get; set; }
        public ulong StartIndex { get; set; }
        public ulong EndIndex { get; set; }
        public List<IFData.Molecule> Molecules { get; set; }
        public IFData.Calibration Calibration { get; set; }
        public IFData.BaselineCorr BaselineCorr { get; set; }

        public double FwhmRange { get; set; }
        public double SearchRange { get; set; }

        public InterpType interpType { get; set; }
        
        #endregion

        public enum InterpType
        {
            Polynomial,
            Spline,
            PCHIP
        }

        #region Methods

        /// <summary>
        /// Loads the contents of an IFD file into the workspace.
        /// </summary>
        /// <param name="path">Path to the IFD file.</param>
        [Obsolete]
        public void LoadIFDFile(string path)
        {
            var rootElement = IFDFile.Open(path);

            SpectralData = IFDFile.ReadRawData(rootElement);
            StartIndex = IFDFile.ReadStartIndex(rootElement);
            EndIndex = IFDFile.ReadEndIndex(rootElement);
            Molecules = IFDFile.ReadMolecules(rootElement);
            Calibration = IFDFile.ReadCalibration(rootElement);
            BaselineCorr = IFDFile.ReadBackgroundCorr(rootElement);
        }


        public void CorrectBaseline()
        {
            //TODO: can be shortened/optimized

            //TODO: check is all necessary data are already loaded in the workspace. throw exception if not
            if (SpectralData == null || BaselineCorr == null)
            {
                throw new WorkspaceNotDefinedException("Required data are not present in workspace."); //TODO: specify more precisely
            }

            int massAxisLength = SpectralData.RawLength;

            //TODO: Evaluating the bg correction for the whole range might be useless. Specifiyng a mass range would make sense.
            //TODO: crop the mass axis to the last specified point of the baseline? might break usefulness.
            PPInterpolation baselineFit = new PPInterpolation(BaselineCorr.XAxis.ToArray(), BaselineCorr.YAxis.ToArray(), PPInterpolation.PPType.PCHIP);    //TODO: write for more interp types

            MathNet.Numerics.LinearAlgebra.Vector<double> baseline = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.DenseOfArray(baselineFit.Evaluate(SpectralData.RawMassAxis.ToArray()));
            MathNet.Numerics.LinearAlgebra.Vector<double> correctedSignal = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(massAxisLength, 0);

            for (int i = 0; i < massAxisLength; i++)
            {
                correctedSignal[i] = SpectralData.RawSignalAxis[i] - baseline[i];
            }

            SpectralData.PureSignalAxis = correctedSignal;
        }

        public void CorrectMassOffset()
        {
            //TODO: this can be cleaned/optimized, it is an ugly mess right now
            //TODO: check if required data are present

            int massAxisLength = SpectralData.RawLength;    //TODO: might want to use the cropped length as well
            int fitDataLength = Calibration.COMList.Count;

            // generate the x-axis for the fit
            double xAxisMin = SpectralData.RawMassAxis.Min();
            double xAxisMax = SpectralData.RawMassAxis.Max();

            int xAxisLength = (int)((xAxisMax - xAxisMin) / 0.01);   //TODO: allow to set granularity of the x axis
            MathNet.Numerics.LinearAlgebra.Vector<double> xAxis = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(xAxisLength);

            for (int i = 0; i < xAxisLength; i++)
            {
                xAxis[i] = xAxisMin + i * 0.01;
            }

            // Will hold evaluated data from the first fit and serve as x-axis (yes, x-axis) in the second fit. Because the fit needs to be reversed.
            double[] yAxis = new double[xAxisLength];

            PPInterpolation massOffset = new PPInterpolation(Calibration.COMList.ToArray(), Calibration.MassOffsetList.ToArray(), PPInterpolation.PPType.PCHIP);

            // Evaluation of the fit at positions given by generated x-axis and storing the correction in YAxis vector.
            yAxis = massOffset.Evaluate(xAxis.ToArray());   //TODO: this is most likely broken because it can not extrapolate right now

            for (int i = 0; i < xAxisLength; i++)
            {
                yAxis[i] += xAxis[i];
            }

            // This is a sanity check for YAxis monotonicity. That is required by the second fit, where the YAxis server as x-axis. Matlab does this without telling.
            // This could be solved in a more safe way, if the user would actually specify endindex of the raw data, so the uncalibrated range can get cut out.
            double yAxisMax = yAxis.Max();
            int yAxisMaxIndex = Array.BinarySearch(yAxis, yAxisMax);
            double[] yAxisNew = new double[yAxisMaxIndex + 1];  //TODO: check if the length correct
            Array.Copy(yAxis, yAxisNew, yAxisNew.Length);

            // Fit to generate corrected mass axis. Note that the X and Y axis are inverted. For this to be correct, it must be corrected in the IFD file generation first.
            PPInterpolation massOffset2 = new PPInterpolation(yAxisNew, xAxis.ToArray(), PPInterpolation.PPType.PCHIP);    //TODO: i think those vectors can have different lengths now

            double[] correctedMassAxis = new double[yAxis.Length];

            Console.WriteLine(SpectralData.RawMassAxis.Last());
            Console.WriteLine(yAxisNew.Last());

            correctedMassAxis = massOffset2.Evaluate(SpectralData.RawMassAxis.ToArray());

            //TODO: this is the old C++ evaluation
            //for (int i = 0; SpectralData.RawMassAxis[i] < yAxisNew.Last(); i++)
            //{
            //    correctedMassAxis[i] = massOffset2.Evaluate(SpectralData.RawMassAxis[i]);
            //}

            // TODO: It might happen due to the monotonicity check, that during the second evaluation we hit a singularity. This cuts off the nonsense data.

            // Store the corrected mass axis in the appropriate place
            SpectralData.MassOffsetCorrAxis = MathNet.Numerics.LinearAlgebra.Double.DenseVector.Build.DenseOfArray(correctedMassAxis);
        }

        public void BuildDesignMatrix()
        {
            designMatrix = new DesignMtrx(SpectralData, Molecules, Calibration);
            designMatrix.Build();
        }

        public void ExtractAbundances()
        {
            LeastSquaresSystem lss = new LeastSquaresSystem(designMatrix.Storage, SpectralData.PureSignalAxis);

            lss = Algorithm.LeaSqrSparseQRHouseholder(lss);

            lss.Solve();
        }

        public void ResolutionFit(InterpType t)
        {
            switch (t)
            {
                case InterpType.Polynomial:
                    int order = 3; // TODO will be defined by user from GUI
                    PolyInterpolation PolyRC = new PolyInterpolation(Calibration.COMList.ToArray(), Calibration.ResolutionList.ToArray(), order);
                    break;
                case InterpType.Spline:
                    PPInterpolation SplineRC = new PPInterpolation(Calibration.COMList.ToArray(), Calibration.ResolutionList.ToArray(), PPInterpolation.PPType.Spline);
                    break;
                case InterpType.PCHIP:
                    PPInterpolation PCHIPRC = new PPInterpolation(Calibration.COMList.ToArray(), Calibration.ResolutionList.ToArray(), PPInterpolation.PPType.PCHIP);
                    break;
                default:
                    throw new Interpolation.InterpolationException("Unknown interpolation type.");
            }

            //TODO: save the object somewhere
        }

        #endregion
    }
}
