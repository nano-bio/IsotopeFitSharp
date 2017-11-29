﻿using System;
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
        public DesignMtrx designMatrix;

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
        public IFData.BaselineCorr BaselineCorrData { get; set; }
        public double[] Abundances { get; private set; }

        public double FwhmRange { get; set; }
        public double SearchRange { get; set; }

        public Interpolation ResolutionInterpolation { get; private set; }
        
        #endregion

        #region Public Methods

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
            BaselineCorrData = IFDFile.ReadBackgroundCorr(rootElement);
        }

        /// <summary>
        /// Calculates baseline corrected signal from raw signal data and baseline correction points. Stores the result in the corresponding Workspace property.
        /// </summary>
        public void CorrectBaseline()
        {
            //TODO: can be shortened/optimized

            if (SpectralData.RawSignalAxis == null) throw new WorkspaceNotDefinedException("Raw signal axis not specified.");
            if (BaselineCorrData.XAxis == null || BaselineCorrData.YAxis == null) throw new WorkspaceNotDefinedException("Baseline correction points not specified.");

            int massAxisLength = SpectralData.RawLength;

            //TODO: Evaluating the bg correction for the whole range might be useless. Specifiyng a mass range would make sense.
            //TODO: crop the mass axis to the last specified point of the baseline? might break usefulness.
            PPInterpolation baselineFit = new PPInterpolation(BaselineCorrData.XAxis.ToArray(), BaselineCorrData.YAxis.ToArray(), PPInterpolation.PPType.PCHIP);    // in the matlab code it is also hard-coded pchip

            MathNet.Numerics.LinearAlgebra.Vector<double> baseline = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.DenseOfArray(baselineFit.Evaluate(SpectralData.RawMassAxis.ToArray()));
            //MathNet.Numerics.LinearAlgebra.Vector<double> correctedSignal = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(massAxisLength, 0);
            double[] correctedSignal = new double[massAxisLength];

            for (int i = 0; i < massAxisLength; i++)
            {
                correctedSignal[i] = SpectralData.RawSignalAxis[i] - baseline[i];
            }

            SpectralData.PureSignalAxis = correctedSignal;
        }

        /// <summary>
        /// Calculates mass axis corrected for mass offset from previously supplied calibration data and stores the result in appropriate location.
        /// </summary>
        /// <param name="interpType">Type of the interpolation to be used.</param>
        /// <param name="order">Order of the polynomial interpolation. This is relevant only for the polynomial interpolation.</param>
        public void CorrectMassOffset(Interpolation.Type interpType, int order)
        {
            //TODO: this can be cleaned/optimized, it is an ugly mess right now

            if (SpectralData.RawMassAxis == null) throw new WorkspaceNotDefinedException("Raw mass axis not specified.");
            if (Calibration.COMList == null || Calibration.MassOffsetList == null) throw new WorkspaceNotDefinedException("Mass offset calibration points not specified.");

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

            // Will hold evaluated data from the first fit and serve as x-axis (yes, x-axis) in the second fit. Because the fit needs to be inverted.
            double[] yAxis = new double[xAxisLength];

            // TODO: temporarily hardcoded for the spline not-a-knot
            //PPInterpolation massOffset = new PPInterpolation(Calibration.COMList.ToArray(), Calibration.MassOffsetList.ToArray(), PPInterpolation.PPType.SplineNotAKnot);

            Interpolation massOffset;

            switch (interpType)
            {
                case Interpolation.Type.Polynomial:
                    massOffset = new PolyInterpolation(Calibration.COMList.ToArray(), Calibration.MassOffsetList.ToArray(), order);
                    break;
                case Interpolation.Type.SplineNatural:
                    throw new NotImplementedException();
                    //break;
                case Interpolation.Type.SplineNotAKnot:
                    massOffset = new PPInterpolation(Calibration.COMList.ToArray(), Calibration.MassOffsetList.ToArray(), PPInterpolation.PPType.SplineNotAKnot);
                    break;
                case Interpolation.Type.PCHIP:
                    massOffset = new PPInterpolation(Calibration.COMList.ToArray(), Calibration.MassOffsetList.ToArray(), PPInterpolation.PPType.PCHIP);
                    break;
                default:
                    throw new Interpolation.InterpolationException("Unknown interpolation type specified.");
            }

            // Evaluation of the fit at positions given by generated x-axis and storing the correction in YAxis vector.
            yAxis = massOffset.Evaluate(xAxis.ToArray());

            for (int i = 0; i < xAxisLength; i++)
            {
                yAxis[i] += xAxis[i];
            }

            // This is a sanity check for YAxis monotonicity. That is required by the second fit, where the YAxis server as x-axis. Matlab does this without telling.
            // This could be solved in a more safe way, if the user would actually specify endindex of the raw data, so the uncalibrated range can get cut out.
            double yAxisMax = yAxis.Max();
            int yAxisMaxIndex = Array.BinarySearch(yAxis, yAxisMax);
            double[] yAxisNew = new double[xAxisLength];  //TODO: check if the length correct  //yAxisMaxIndex + 1
            Array.Copy(yAxis, yAxisNew, yAxisNew.Length);

            // Fit to generate corrected mass axis. Note that the X and Y axis are inverted. For this to be correct, it must be corrected in the IFD file generation first.
            PPInterpolation massOffset2 = new PPInterpolation(yAxisNew, xAxis.ToArray(), PPInterpolation.PPType.PCHIP);    //TODO: i think those vectors can have different lengths now

            double[] correctedMassAxis = new double[yAxis.Length];

            correctedMassAxis = massOffset2.Evaluate(SpectralData.RawMassAxis.ToArray());

            //TODO: this is the old C++ evaluation
            //for (int i = 0; SpectralData.RawMassAxis[i] < yAxisNew.Last(); i++)
            //{
            //    correctedMassAxis[i] = massOffset2.Evaluate(SpectralData.RawMassAxis[i]);
            //}

            // TODO: It might happen due to the monotonicity check, that during the second evaluation we hit a singularity. This cuts off the nonsense data.

            // Store the corrected mass axis in the appropriate place
            SpectralData.MassOffsetCorrAxis = correctedMassAxis;  //MathNet.Numerics.LinearAlgebra.Double.DenseVector.Build.DenseOfArray(correctedMassAxis);
        }

        /// <summary>
        /// Fits the previously supplied resolution calibration data and stores the calibration results.
        /// </summary>
        /// <param name="t">Type of the interpolation ti use.</param>
        /// <param name="order">Order of the polynomial interpolation. This is relevant only for the polynomial interpolation.</param>
        public void ResolutionFit(Interpolation.Type t, int order)
        {
            switch (t)
            {
                case Interpolation.Type.Polynomial:
                    //int order = 3; // TODO will be defined by user from GUI
                    ResolutionInterpolation = new PolyInterpolation(Calibration.COMList.ToArray(), Calibration.ResolutionList.ToArray(), order); //PolyInterpolation PolyRC
                    break;
                case Interpolation.Type.SplineNatural:
                    throw new NotImplementedException("This has not yet been implemented.");
                //break;
                case Interpolation.Type.SplineNotAKnot:
                    ResolutionInterpolation = new PPInterpolation(Calibration.COMList.ToArray(), Calibration.ResolutionList.ToArray(), PPInterpolation.PPType.SplineNotAKnot); //PPInterpolation SplineRC
                    break;
                case Interpolation.Type.PCHIP:
                    ResolutionInterpolation = new PPInterpolation(Calibration.COMList.ToArray(), Calibration.ResolutionList.ToArray(), PPInterpolation.PPType.PCHIP); //PPInterpolation PCHIPRC
                    break;
                default:
                    throw new Interpolation.InterpolationException("Unknown interpolation type.");
            }
        }

        public void BuildDesignMatrix()
        {
            designMatrix = new DesignMtrx(SpectralData, Molecules, Calibration, ResolutionInterpolation);
            designMatrix.Build();
        }

        public void ExtractAbundances()
        {
            LeastSquaresSystem lss = new LeastSquaresSystem(designMatrix.Storage, null);    // we dont give the observation vector here, because it was already added to the matrix during the build
            
            //TODO: this is ugly, hide this inside the LeastSquaresSystem class
            lss = Algorithm.LeaSqrSparseQRHouseholder(lss);

            lss.Solve();

            Abundances = lss.Solution.ToArray();
        }

        #endregion
    }
}
