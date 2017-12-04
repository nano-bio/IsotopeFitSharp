using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

//using IsotopeFit.Numerics;

namespace IsotopeFit
{
    //[ComVisible(true)]
	/// <summary>
	/// This class is the main entry point for outside usage.
	/// </summary>
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
            MathNet.Numerics.Control.UseNativeMKL();

            SpectralData = new IFData.Spectrum();
            Clusters = new OrderedDictionary();
            Calibration = new IFData.Calibration();
            BaselineCorrData = new IFData.BaselineCorr();
        }

        /// <summary>
        /// Create an IsotopeFit workspace and load an IFD/IFJ file into it.
        /// </summary>
        /// <param name="IFDfile">Path to the IFD file to be loaded.</param>
        public Workspace(string path)
        {
            MathNet.Numerics.Control.UseNativeMKL();

            LoadIFDFile(path);            
        }

        #endregion

        #region Destructors

        #endregion

        #region Properties

        /// <summary>
        /// Object containing spectral data, both raw and calibrated.
        /// </summary>
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
        public OrderedDictionary Abundances { get; private set; }
        //public double[] Abundances { get; private set; }

        public double FwhmRange { get; set; }
        public double SearchRange { get; set; }

        public Interpolation ResolutionInterpolation { get; private set; }
        public DesignMtrx DesignMatrix { get; private set; }

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
            //TODO: make more robust, recognize the filetype

            var rootElement = IFDFile.Open(path);

            SpectralData = IFDFile.ReadRawData(rootElement);
            StartIndex = IFDFile.ReadStartIndex(rootElement);
            EndIndex = IFDFile.ReadEndIndex(rootElement);
            Clusters = IFDFile.ReadMolecules(rootElement);
            Calibration = IFDFile.ReadCalibration(rootElement);
            BaselineCorrData = IFDFile.ReadBackgroundCorr(rootElement);

            IFDLoaded = true;
        }

        /// <summary>
        /// Calculates baseline corrected signal from raw signal data and baseline correction points. Stores the result in the <see cref="Workspace.SpectralData"/>.SignalAxis property.
        /// </summary>
        /// <remarks>
        /// This method uses the PCHIP interpolation to obtain the baseline values, which are stored in the <see cref="Workspace.SpectralData"/>.Baseline property.
        /// Optional arguments are meant to be used by GUI to ease the process of loading data to the <see cref="Workspace"/>. If not supplied, the funcion will use previously stored values.
        /// </remarks>
        /// <param name="xAxis">Optional x-axis for the baseline correction.</param>
        /// <param name="yAxis">Optional y-axis for the baseline correction.</param>
        /// <exception cref="WorkspaceNotDefinedException">Thrown when some of the required data have not been loaded in the <see cref="Workspace"/>.</exception>
        public void CorrectBaseline(double[] xAxis = null, double[] yAxis = null)
        {
            // raw spectral data are required
            if (SpectralData.RawSignalAxis == null || SpectralData.RawMassAxis == null) throw new WorkspaceNotDefinedException("Raw spectral data not specified.");
            if (SpectralData.RawSignalAxis.Length != SpectralData.RawMassAxis.Length) throw new WorkspaceNotDefinedException("Supplied spectral data axis have different lengths.");

            // if there are optional background correction points specified, save them to the workspace and use those
            if (xAxis != null) BaselineCorrData.XAxis = xAxis;
            if (yAxis != null) BaselineCorrData.YAxis = yAxis;

            // safety check, they can still be null in come use cases
            if (BaselineCorrData.XAxis == null || BaselineCorrData.YAxis == null) throw new WorkspaceNotDefinedException("Baseline correction points not specified.");
            if (BaselineCorrData.XAxis.Length != BaselineCorrData.YAxis.Length) throw new WorkspaceNotDefinedException("Supplied baseline correction point arrays have different lengths.");

            int massAxisLength = SpectralData.RawLength;

            //TODO: Evaluating the bg correction for the whole range might be useless. Specifiyng a mass range would make sense.
            //TODO: crop the mass axis to the last specified point of the baseline? might break usefulness.
            PPInterpolation baselineFit = new PPInterpolation(BaselineCorrData.XAxis, BaselineCorrData.YAxis, PPInterpolation.PPType.PCHIP);    // in the matlab code it is also hard-coded pchip

            SpectralData.Baseline = new double[massAxisLength];
            SpectralData.SignalAxis = new double[massAxisLength];

            for (int i = 0; i < massAxisLength; i++)
            {
                SpectralData.Baseline[i] = baselineFit.Evaluate(SpectralData.RawMassAxis[i]);
                SpectralData.SignalAxis[i] = SpectralData.RawSignalAxis[i] - SpectralData.Baseline[i];
            }
        }

        /// <summary>
        /// Calculates mass axis corrected for mass offset from previously supplied calibration data and stores the result in the <see cref="Workspace.SpectralData"/>.MassAxis property.
        /// </summary>
        /// <remarks>
        /// In case of polynomial inerpolation, the order is considered the highest exponent value of the desired interpolating polynomial, e.g. order = 2 will produce quadratic polynomial.
        /// Optional arguments are meant to be used by GUI to ease the process of loading data to the <see cref="Workspace"/>. If not supplied, the method will use previously stored values.
        /// </remarks>
        /// <param name="interpType">Type of the interpolation to be used.</param>
        /// <param name="order">Order of the polynomial interpolation. Ignored if any other interpolation type is selected.</param>
        /// <param name="comList">Optional array of centre-of-mass points to be used for the correction.</param>
        /// <param name="massOffsetList">Optional array of mass offset points to be used for the correction.</param>
        /// <exception cref="WorkspaceNotDefinedException">Thrown when some of the required data have not been loaded in the <see cref="Workspace"/> or their lengths are different.</exception>
        public void CorrectMassOffset(Interpolation.Type interpType, int order, double[] comList = null, double[] massOffsetList = null)
        {
            //TODO: this can be cleaned/optimized, it is an ugly mess right now

            // required data
            if (SpectralData.RawMassAxis == null) throw new WorkspaceNotDefinedException("Raw mass axis not specified.");

            // there are three use cases
            if (!IFDLoaded && comList == null && massOffsetList == null)    // this is the IFJFile or GUI use case, when the COMList and MassOffsetList must be built
            {
                if (Calibration.NameList.Count == 0) throw new WorkspaceNotDefinedException("Calibration namelist empty.");

                int n = Calibration.NameList.Count;

                Calibration.COMList = new double[n];
                Calibration.MassOffsetList = new double[n];

                for (int i = 0; i < n; i++)
                {
                    Calibration.COMList[i] = (Clusters[Calibration.NameList[i]] as IFData.Cluster).CentreOfMass;
                    Calibration.MassOffsetList[i] = (Clusters[Calibration.NameList[i]] as IFData.Cluster).MassOffset;
                }
            }
            else    // this is the second use case for GUI, because only the GUI will supply the optional arguments if needed. IFDFile case just passes over those ifs
            {
                if (comList != null) Calibration.COMList = comList;
                if (massOffsetList != null) Calibration.MassOffsetList = massOffsetList;
            }

            // but we still need to check, because in the third use case (IFDFile) the data must still be valid
            if (Calibration.COMList == null || Calibration.MassOffsetList == null) throw new WorkspaceNotDefinedException("Mass offset calibration points not specified.");
            if (Calibration.COMList.Length != Calibration.MassOffsetList.Length) throw new WorkspaceNotDefinedException("Supplied mass offset correction point arrays have different lengths.");

            //int massAxisLength = SpectralData.RawLength;    //TODO: might want to use the cropped length as well
            //int fitDataLength = Calibration.COMList.Length;

            // generate the x-axis for the fit
            double xAxisMin = SpectralData.RawMassAxis.Min();
            double xAxisMax = SpectralData.RawMassAxis.Max();

            int xAxisLength = (int)((xAxisMax - xAxisMin) / 0.01);   //TODO: allow to set granularity of the x axis?
            //MathNet.Numerics.LinearAlgebra.Vector<double> xAxis = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(xAxisLength);
            double[] xAxis = new double[xAxisLength];

            for (int i = 0; i < xAxisLength; i++)
            {
                xAxis[i] = xAxisMin + i * 0.01;
            }

            // first interpolation
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

            // Will hold evaluated data from the first fit and serve as x-axis (yes, x-axis) in the second fit. Because the fit needs to be inverted.
            // Evaluation of the fit at positions given by generated x-axis and storing the correction in YAxis vector.
            double[] yAxis = new double[xAxisLength];
            //yAxis = massOffset.Evaluate(xAxis.ToArray());

            for (int i = 0; i < xAxisLength; i++)
            {
                //yAxis[i] += xAxis[i];
                yAxis[i] = massOffset.Evaluate(xAxis[i]) + xAxis[i];
            }

            // This is a sanity check for YAxis monotonicity. That is required by the second fit, where the YAxis server as x-axis. Matlab does this without telling.
            // This could be solved in a more safe way, if the user would actually specify endindex of the raw data, so the uncalibrated range can get cut out.
            double yAxisMax = yAxis.Max();  //TODO: this is not yet robust
            int yAxisMaxIndex = Array.BinarySearch(yAxis, yAxisMax);
            double[] yAxisNew = new double[xAxisLength];  //TODO: check if the length correct  //yAxisMaxIndex + 1
            Array.Copy(yAxis, yAxisNew, yAxisNew.Length);

            // Fit to generate corrected mass axis. Note that the X and Y axis are inverted. For this to be correct, it must be corrected in the IFD file generation first.
            PPInterpolation massOffset2 = new PPInterpolation(yAxisNew, xAxis.ToArray(), PPInterpolation.PPType.PCHIP);    //TODO: i think those vectors can have different lengths now

            double[] correctedMassAxis = new double[yAxis.Length];

            correctedMassAxis = massOffset2.Evaluate(SpectralData.RawMassAxis.ToArray());

            // TODO: It might happen due to the monotonicity check, that during the second evaluation we hit a singularity. This cuts off the nonsense data.

            // Store the corrected mass axis in the appropriate place
            SpectralData.MassAxis = correctedMassAxis;
        }

        /// <summary>
        /// Fits the previously supplied resolution calibration data and stores the calibration results in the <see cref="Workspace.ResolutionInterpolation"/> property.
        /// </summary>
        /// <remarks>
        /// In case of polynomial inerpolation, the order is considered the highest exponent value of the desired interpolating polynomial, e.g. order = 2 will produce quadratic polynomial.
        /// Optional arguments are meant to be used by GUI to ease the process of loading data to the <see cref="Workspace"/>. If not supplied, the method will use previously stored values.
        /// </remarks>
        /// <param name="t">Type of the interpolation to use.</param>
        /// <param name="order">Order of the polynomial interpolation. This is relevant only for the polynomial interpolation.</param>
        /// <param name="comList">Optional array of centre-of-mass points to be used for the calibration.</param>
        /// <param name="resolutionList">Optional array of resolution points to be used for the calibration.</param>
        /// <exception cref="WorkspaceNotDefinedException">Thrown when some of the required data have not been loaded in the <see cref="Workspace"/> or their lengths are different.</exception>
        public void ResolutionFit(Interpolation.Type t, int order, double[] comList = null, double[] resolutionList = null)
        {
            // there are three use cases
            if (!IFDLoaded && comList == null && resolutionList == null) // this covers the use cases for IFJFile and GUI, when the COMList and ResolutionList have to be built.
            {
                if (Calibration.NameList.Count == 0) throw new WorkspaceNotDefinedException("Calibration namelist empty.");

                int n = Calibration.NameList.Count;

                Calibration.COMList = new double[n];
                Calibration.ResolutionList = new double[n];

                for (int i = 0; i < n; i++)
                {
                    Calibration.COMList[i] = (Clusters[Calibration.NameList[i]] as IFData.Cluster).CentreOfMass;
                    Calibration.ResolutionList[i] = (Clusters[Calibration.NameList[i]] as IFData.Cluster).Resolution;
                }
            }
            else    // this is the second use case for GUI, because only the GUI will supply the optional arguments if needed. IFDFile case just passes over those ifs
            {
                if (comList != null) Calibration.COMList = comList;
                if (resolutionList != null) Calibration.ResolutionList = resolutionList;
            }

            // safety check, because in the third use case (IFDFile) the data must still be valid
            if (Calibration.COMList == null || Calibration.ResolutionList == null) throw new WorkspaceNotDefinedException("Resolution calibration points not specified.");
            if (Calibration.COMList.Length != Calibration.ResolutionList.Length) throw new WorkspaceNotDefinedException("Supplied resolution calibration point arrays have different lengths.");

            switch (t)
            {
                case Interpolation.Type.Polynomial:
                    ResolutionInterpolation = new PolyInterpolation(Calibration.COMList, Calibration.ResolutionList, order);
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
        /// Builds the design matrix from currently supplied data and stores the result in the <see cref="Workspace.DesignMatrix"/> property.
        /// </summary>
        public void BuildDesignMatrix()
        {
            // TODO: merge this with the fit abundances maybe? If we do the design matrix updating, we will still need this function.
            DesignMatrix = new DesignMtrx(SpectralData, Clusters, Calibration, ResolutionInterpolation);
            DesignMatrix.Build();
        }

        /// <summary>
        /// Performs the fit of the data and stores the result in the <see cref="Workspace.Abundances"/> property.
        /// </summary>
        public void FitAbundances()
        {
            LeastSquaresSystem lss = new LeastSquaresSystem(DesignMatrix.Storage, null);    // we dont give the observation vector here, because it was already added to the matrix during the build
            
            //TODO: this is ugly, hide this inside the LeastSquaresSystem class
            lss = Algorithm.LeaSqrSparseQRHouseholder(lss);

            lss.Solve();

            Abundances = new OrderedDictionary(Clusters.Count);

            for (int i = 0; i < Clusters.Count; i++)
            {
                Abundances.Add((Clusters[i] as IFData.Cluster).Name, lss.Solution[i]);
            }
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
