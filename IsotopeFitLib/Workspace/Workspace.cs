using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

using IsotopeFit.Numerics;

namespace IsotopeFit
{
    //[ComVisible(true)]
    public partial class Workspace
    {
        #region Fields
        

        #endregion

        #region Constructors

        /// <summary>
        /// Create an empty IsotopeFit workspace.
        /// </summary>
        public Workspace()
        {
            SpectralData = new IFData.Spectrum();
            Clusters = new OrderedDictionary();
            Calibration = new IFData.Calibration();
            BaselineCorrData = new IFData.BaselineCorr();
        }

        /// <summary>
        /// Create an IsotopeFit workspace and load an IFD file into it.
        /// </summary>
        /// <param name="IFDfile">Path to the IFD file to be loaded.</param>
        public Workspace(string path)
        {
            LoadIFDFile(path);
            IFDLoaded = true;
        }

        #endregion

        #region Destructors

        #endregion

        #region Properties

        public IFData.Spectrum SpectralData { get; set; }

        /// <summary>
        /// Crop start index. Not implemented at the moment.
        /// </summary>
        public int StartIndex { get; set; }

        /// <summary>
        /// Crop end index. Not implemented at the moment.
        /// </summary>
        public int EndIndex { get; set; }

        /// <summary>
        /// Object containing the clusters to be fitted.
        /// </summary>
        public OrderedDictionary Clusters { get; set; } //TODO: maybe I should add a convenience functions for adding a cluster to the dictionary, that way I can check the inputs immediately

        /// <summary>
        /// Object containing the data necessary for mass offset correction, resolution fit and peak shape.
        /// </summary>
        public IFData.Calibration Calibration { get; set; }

        /// <summary>
        /// Object containing the data necessary for baseline correction.
        /// </summary> 
        public IFData.BaselineCorr BaselineCorrData { get; set; }

        /// <summary>
        /// Object containing the fitted cluster abundances.
        /// </summary>
        public double[] Abundances { get; private set; }

        public double FwhmRange { get; set; }
        public double SearchRange { get; set; }

        public Interpolation ResolutionInterpolation { get; private set; }
        public DesignMtrx DesignMatrix { get; internal set; }

        public WorkspaceStatus Status { get; private set; }

        internal bool IFDLoaded { get; private set; }

        #endregion

        #region Public Methods

        /// <summary>
        /// Loads the contents of an IFD file into the workspace.
        /// </summary>
        /// <param name="path">Path to the IFD file.</param>
        public void LoadIFDFile(string path)
        {
            var rootElement = IFDFile.Open(path);

            SpectralData = IFDFile.ReadRawData(rootElement);
            StartIndex = IFDFile.ReadStartIndex(rootElement);
            EndIndex = IFDFile.ReadEndIndex(rootElement);
            Clusters = IFDFile.ReadMolecules(rootElement);
            Calibration = IFDFile.ReadCalibration(rootElement);
            BaselineCorrData = IFDFile.ReadBackgroundCorr(rootElement);
        }

        /// <summary>
        /// Calculates baseline corrected signal from raw signal data and baseline correction points. Stores the result in the Workspace.SpectralData.SignalAxis property.
        /// </summary>
        /// <remarks>Uses the PCHIP interpolation.</remarks>
        public void CorrectBaseline()
        {
            //TODO: can be shortened/optimized

            if (SpectralData.RawSignalAxis == null || SpectralData.RawMassAxis == null) throw new WorkspaceNotDefinedException("Raw spectral data not specified.");
            if (SpectralData.RawSignalAxis.Length != SpectralData.RawMassAxis.Length) throw new WorkspaceNotDefinedException("Supplied spectral data axis have different lengths.");

            if (BaselineCorrData.XAxis == null || BaselineCorrData.YAxis == null) throw new WorkspaceNotDefinedException("Baseline correction points not specified.");
            if (BaselineCorrData.XAxis.Length != BaselineCorrData.YAxis.Length) throw new WorkspaceNotDefinedException("Supplied baseline correction point arrays have different lengths.");

            int massAxisLength = SpectralData.RawLength;

            //TODO: Evaluating the bg correction for the whole range might be useless. Specifiyng a mass range would make sense.
            //TODO: crop the mass axis to the last specified point of the baseline? might break usefulness.
            PPInterpolation baselineFit = new PPInterpolation(BaselineCorrData.XAxis.ToArray(), BaselineCorrData.YAxis.ToArray(), PPInterpolation.PPType.PCHIP);    // in the matlab code it is also hard-coded pchip
                        
            SpectralData.SignalAxis = new double[massAxisLength];

            for (int i = 0; i < massAxisLength; i++)
            {
                SpectralData.SignalAxis[i] = SpectralData.RawSignalAxis[i] - baselineFit.Evaluate(SpectralData.RawMassAxis[i]);
            }
        }

        /// <summary>
        /// Calculates mass axis corrected for mass offset from previously supplied calibration data and stores the result in the Workspace.SpectralData.MassAxis property.
        /// </summary>
        /// <param name="interpType">Type of the interpolation to be used.</param>
        /// <param name="order">Order of the polynomial interpolation. This is relevant only for the polynomial interpolation.</param>
        public void CorrectMassOffset(Interpolation.Type interpType, int order)
        {
            //TODO: this can be cleaned/optimized, it is an ugly mess right now
            //TODO: implement check for the Calibration.Namelist and then decide which strategy will we follow

            if (SpectralData.RawMassAxis == null) throw new WorkspaceNotDefinedException("Raw mass axis not specified.");

            // check if we have loaded an IFD file, then we can continue. Otherwise we are working with IFJ or GUI and we need to create the lists.
            if (!IFDLoaded)
            {
                if (Calibration.NameList.Count == 0) throw new WorkspaceNotDefinedException("Calibration namelist empty.");

                int n = Calibration.NameList.Count;

                double[] comList = new double[n];
                double[] massOffsetList = new double[n];

                for (int i = 0; i < n; i++)
                {
                    comList[i] = (Clusters[Calibration.NameList[i]] as IFData.Cluster).CentreOfMass;
                    massOffsetList[i] = (Clusters[Calibration.NameList[i]] as IFData.Cluster).MassOffset;
                }

                Calibration.COMList = comList;
                Calibration.MassOffsetList = massOffsetList;
            }

            if (Calibration.COMList == null || Calibration.MassOffsetList == null) throw new WorkspaceNotDefinedException("Mass offset calibration points not specified.");

            int massAxisLength = SpectralData.RawLength;    //TODO: might want to use the cropped length as well
            int fitDataLength = Calibration.COMList.Length;

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
            SpectralData.MassAxis = correctedMassAxis;  //MathNet.Numerics.LinearAlgebra.Double.DenseVector.Build.DenseOfArray(correctedMassAxis);
        }

        /// <summary>
        /// Fits the previously supplied resolution calibration data and stores the calibration results in the Workspace.ResolutionInterpolation property.
        /// </summary>
        /// <param name="t">Type of the interpolation to use.</param>
        /// <param name="order">Order of the polynomial interpolation. This is relevant only for the polynomial interpolation.</param>
        public void ResolutionFit(Interpolation.Type t, int order)
        {
            // check if we have loaded an IFD file, then we can continue. Otherwise we are working with IFJ or GUI and we need to create the lists.
            if (!IFDLoaded)
            {
                if (Calibration.NameList.Count == 0) throw new WorkspaceNotDefinedException("Calibration namelist empty.");

                int n = Calibration.NameList.Count;

                double[] comList = new double[n];
                double[] massOffsetList = new double[n];

                for (int i = 0; i < n; i++)
                {
                    comList[i] = (Clusters[Calibration.NameList[i]] as IFData.Cluster).CentreOfMass;
                    massOffsetList[i] = (Clusters[Calibration.NameList[i]] as IFData.Cluster).MassOffset;
                }

                Calibration.COMList = comList;
                Calibration.MassOffsetList = massOffsetList;
            }

            if (Calibration.COMList == null || Calibration.MassOffsetList == null) throw new WorkspaceNotDefinedException("Resolution calibration points not specified.");

            switch (t)
            {
                case Interpolation.Type.Polynomial:
                    ResolutionInterpolation = new PolyInterpolation(Calibration.COMList.ToArray(), Calibration.ResolutionList.ToArray(), order);
                    break;
                case Interpolation.Type.SplineNatural:
                    throw new NotImplementedException("Natural spline interpolation has not yet been implemented.");
                //break;
                case Interpolation.Type.SplineNotAKnot:
                    ResolutionInterpolation = new PPInterpolation(Calibration.COMList.ToArray(), Calibration.ResolutionList.ToArray(), PPInterpolation.PPType.SplineNotAKnot);
                    break;
                case Interpolation.Type.PCHIP:
                    ResolutionInterpolation = new PPInterpolation(Calibration.COMList.ToArray(), Calibration.ResolutionList.ToArray(), PPInterpolation.PPType.PCHIP);
                    break;
                default:
                    throw new Interpolation.InterpolationException("Unknown interpolation type.");
            }
        }

        /// <summary>
        /// Builds the design matrix from currently supplied data and stores the result in the Workspace.DesignMatrix property.
        /// </summary>
        public void BuildDesignMatrix()
        {
            // TODO: merge this with the fit abundances maybe? If we do the design matrix updating, we will still need this function.
            DesignMatrix = new DesignMtrx(SpectralData, Clusters, Calibration, ResolutionInterpolation);
            DesignMatrix.Build();
        }

        /// <summary>
        /// Performs the fit of the data and stores the result in the Workspace.Abundances property.
        /// </summary>
        public void FitAbundances()
        {
            LeastSquaresSystem lss = new LeastSquaresSystem(DesignMatrix.Storage, null);    // we dont give the observation vector here, because it was already added to the matrix during the build
            
            //TODO: this is ugly, hide this inside the LeastSquaresSystem class
            lss = Algorithm.LeaSqrSparseQRHouseholder(lss);

            lss.Solve();

            Abundances = lss.Solution.ToArray();
        }

        #endregion

        /// <summary>
        /// Class for storing the Workspace status flags and message log for the user.
        /// </summary>
        public class WorkspaceStatus
        {
            internal WorkspaceStatus() { }

            public bool BaselineCorrected { get; internal set; }
            public bool MassOffsetCorrected { get; internal set; }
            public bool DesignMatrixBuilt { get; internal set; }
            public bool AbundancesFitted { get; internal set; }

            public List<double> MessageLog { get; internal set; }
        }
    }
}
