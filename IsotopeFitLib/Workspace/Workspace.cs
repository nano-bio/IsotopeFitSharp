using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

namespace IsotopeFit
{
	/// <summary>
	/// This class is the main interface for outside usage.
	/// </summary>
    public partial class Workspace
    {
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

        #region Properties

        /// <summary>
        /// Object containing spectral data, both raw and calibrated.
        /// </summary>
        public IFData.Spectrum SpectralData { get; set; }

        /// <summary>
        /// Object containing the clusters, abundance of which is to be calculated. See also <see cref="IFData.Cluster"/>.
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

        ///// <summary>
        ///// Object containing the fitted cluster abundances.
        ///// </summary>
        //public OrderedDictionary Abundances { get; private set; }
        //public double[] Abundances { get; private set; }
        
        /// <summary>
        /// Object containing the calculated design matrix for current cluster system.
        /// </summary>
        public DesignMtrx DesignMatrix { get; private set; }

        //public WorkspaceStatus Status { get; private set; }

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
            //SpectralData.CropStartIndex = IFDFile.ReadStartIndex(rootElement);
            //SpectralData.CropEndIndex = IFDFile.ReadEndIndex(rootElement);
            Clusters = IFDFile.ReadMolecules(rootElement);
            Calibration = IFDFile.ReadCalibration(rootElement);
            BaselineCorrData = IFDFile.ReadBackgroundCorr(rootElement);

            IFDLoaded = true;
        }

        /// <summary>
        /// Calculates baseline corrected signal from raw signal data and baseline correction points. Stores the result in the <see cref="Workspace.SpectralData"/>.SignalAxis property.
        /// </summary>
        /// <remarks>
        /// <para>This method uses the PCHIP interpolation to obtain the baseline values, which are stored in the <see cref="Workspace.SpectralData"/>.Baseline property.</para>
        /// <para>Optional arguments <paramref name="xAxis"/> and <paramref name="yAxis"/> are meant to ease the process of loading data to the <see cref="Workspace"/>.
        /// If not supplied, the funcion will use previously stored values.</para>
        /// </remarks>
        /// <param name="xAxis">Optional x-axis for the baseline correction.</param>
        /// <param name="yAxis">Optional y-axis for the baseline correction.</param>
        /// <exception cref="WorkspaceException">Thrown when some of the required data have not been loaded in the <see cref="Workspace"/>.</exception>
        public void CorrectBaseline(double[] xAxis = null, double[] yAxis = null)
        {
            // raw spectral data are required
            if (SpectralData.RawSignalAxis == null || SpectralData.RawMassAxis == null) throw new WorkspaceException("Raw spectral data not specified.");
            if (SpectralData.RawSignalAxis.Length != SpectralData.RawMassAxis.Length) throw new WorkspaceException("Supplied spectral data axis have different lengths.");

            // if there are optional background correction points specified, save them to the workspace and use those
            if (xAxis != null) BaselineCorrData.XAxis = xAxis;
            if (yAxis != null) BaselineCorrData.YAxis = yAxis;

            // safety check, they can still be null in come use cases
            if (BaselineCorrData.XAxis == null || BaselineCorrData.YAxis == null) throw new WorkspaceException("Baseline correction points not specified.");
            if (BaselineCorrData.XAxis.Length != BaselineCorrData.YAxis.Length) throw new WorkspaceException("Supplied baseline correction point arrays have different lengths.");

            int massAxisLength = SpectralData.RawLength;

            //TODO: Evaluating the bg correction for the whole range might be useless. Specifiyng a mass range would make sense.
            //TODO: crop the mass axis to the last specified point of the baseline? might break usefulness.
            //PPInterpolation baselineFit = new PPInterpolation(BaselineCorrData.XAxis, BaselineCorrData.YAxis, PPInterpolation.PPType.PCHIP);    // in the matlab code it is also hard-coded pchip
            BaselineCorrData.BaselineInterpolation = new PPInterpolation(BaselineCorrData.XAxis, BaselineCorrData.YAxis, PPInterpolation.PPType.PCHIP);    // in the matlab code it is also hard-coded pchip

            //SpectralData.Baseline = new double[massAxisLength];
            SpectralData.SignalAxis = new double[massAxisLength];

            for (int i = 0; i < massAxisLength; i++)
            {
                //SpectralData.Baseline[i] = BaselineCorrData.BaselineInterpolation.Evaluate(SpectralData.RawMassAxis[i]);
                //SpectralData.SignalAxis[i] = SpectralData.RawSignalAxis[i] - SpectralData.Baseline[i];
                SpectralData.SignalAxis[i] = SpectralData.RawSignalAxis[i] - BaselineCorrData.BaselineInterpolation.Evaluate(SpectralData.RawMassAxis[i]);
            }
        }

        ///// <summary>
        ///// Function to set the mass axis cropping.
        ///// </summary>
        ///// <remarks>
        ///// <para>The spectral data outside the specified interval will not be discarded, merely ignored in later calculations.</para>
        ///// </remarks>
        ///// <param name="startIndex">Start index of the interval to select.</param>
        ///// <param name="endIndex">End index of the interval to select.</param>
        //public void CropMassAxis(int startIndex, int endIndex)
        //{
        //    //TODO: when is the best time to crop? do we know even before the baseline correction what interval do we want?

        //    if (startIndex >= endIndex) throw new WorkspaceException("Crop start index must be less than the end index.");
        //    if (startIndex < 0 || startIndex <= SpectralData.RawLength) throw new WorkspaceException("Crop start index must be non-negative and less than raw mass axis length.");
        //    if (endIndex < 0 || endIndex <= SpectralData.RawLength) throw new WorkspaceException("Crop start index must be non-negative and less than raw mass axis length.");

        //    SpectralData.CropStartIndex = startIndex;
        //    SpectralData.CropEndIndex = endIndex;

        //    SpectralData.MassAxisCrop = new double[endIndex - startIndex + 1];
        //    SpectralData.SignalAxisCrop = new double[endIndex - startIndex + 1];

        //    SpectralData.Cropped = true;
        //    Array.Copy(SpectralData.RawMassAxis, startIndex, SpectralData.MassAxisCrop, 0, endIndex - startIndex + 1);
        //    Array.Copy(SpectralData.RawSignalAxis, startIndex, SpectralData.SignalAxisCrop, 0, endIndex - startIndex + 1);
        //}

        /// <summary>
        /// Calculates mass axis corrected for mass offset from previously supplied calibration data and stores the result in the <see cref="Workspace.SpectralData"/>.MassAxis property.
        /// </summary>
        /// <remarks>
        /// <para>In case of polynomial inerpolation, the order is considered the highest exponent value of the desired interpolating polynomial, e.g. order = 2 will produce quadratic polynomial.</para>
        /// <para>Optional arguments <paramref name="comList"/> and <paramref name="massOffsetList"/> are meant to ease the process of loading data to the <see cref="Workspace"/>.
        /// If not supplied, the method will use previously stored values.</para>
        /// </remarks>
        /// <param name="interpType">Type of the interpolation to be used.</param>
        /// <param name="order">Order of the polynomial interpolation. Is marked as optional, but is required if polynomial interpolation is selected.</param>
        /// <param name="comList">Optional array of centre-of-mass points to be used for the correction.</param>
        /// <param name="massOffsetList">Optional array of mass offset points to be used for the correction.</param>
        /// <exception cref="WorkspaceException">Thrown when some of the required data have not been loaded in the <see cref="Workspace"/> or their lengths are different.</exception>
        public void CorrectMassOffset(Interpolation.Type interpType, int order = -1, double[] comList = null, double[] massOffsetList = null)
        {
            //TODO: this can be cleaned/optimized, it is an ugly mess right now
            //TODO: add autocrop after the last calibrated peak and maybe before the first one as well

            // if a polynomial fir is requested, check that the order was supplied as well
            if (interpType == Interpolation.Type.Polynomial && order < 0) throw new WorkspaceException("The polynomial order must be specified for the polynomial interpolation type.");

            // required data
            if (SpectralData.RawMassAxis == null) throw new WorkspaceException("Raw mass axis not specified.");

            // there are three use cases
            if (!IFDLoaded && comList == null && massOffsetList == null)    // this is the IFJFile or GUI use case, when the COMList and MassOffsetList must be built
            {
                if (Calibration.NameList.Count == 0) throw new WorkspaceException("Calibration namelist empty.");

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
            if (Calibration.COMList == null || Calibration.MassOffsetList == null) throw new WorkspaceException("Mass offset calibration points not specified.");
            if (Calibration.COMList.Length != Calibration.MassOffsetList.Length) throw new WorkspaceException("Supplied mass offset correction point arrays have different lengths.");

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
            //Interpolation massOffset;

            switch (interpType)
            {
                case Interpolation.Type.Polynomial:
                    Calibration.MassOffsetInterp = new PolyInterpolation(Calibration.COMList, Calibration.MassOffsetList, order);
                    break;
                case Interpolation.Type.SplineNatural:
                    throw new NotImplementedException();
                    //break;
                case Interpolation.Type.SplineNotAKnot:
                    Calibration.MassOffsetInterp = new PPInterpolation(Calibration.COMList, Calibration.MassOffsetList, PPInterpolation.PPType.SplineNotAKnot);
                    break;
                case Interpolation.Type.PCHIP:
                    Calibration.MassOffsetInterp = new PPInterpolation(Calibration.COMList, Calibration.MassOffsetList, PPInterpolation.PPType.PCHIP);
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
                yAxis[i] = Calibration.MassOffsetInterp.Evaluate(xAxis[i]) + xAxis[i];
            }

            // This is a sanity check for YAxis monotonicity. That is required by the second fit, where the YAxis server as x-axis. Matlab does this without telling.
            // This could be solved in a more safe way, if the user would actually specify endindex of the raw data, so the uncalibrated range can get cut out.
            double yAxisMax = yAxis.Max();  //TODO: this is not yet robust
            int yAxisMaxIndex = Array.BinarySearch(yAxis, yAxisMax);
            double[] yAxisNew = new double[xAxisLength];  //TODO: check if the length correct  //yAxisMaxIndex + 1
            Array.Copy(yAxis, yAxisNew, yAxisNew.Length);

            // Fit to generate corrected mass axis. Note that the X and Y axis are inverted. For this to be correct, it must be corrected in the IFD file generation first.
            PPInterpolation correctMassInterp = new PPInterpolation(yAxisNew, xAxis, PPInterpolation.PPType.PCHIP);    //TODO: i think those vectors can have different lengths now

            double[] correctedMassAxis = new double[yAxis.Length];

            correctedMassAxis = correctMassInterp.Evaluate(SpectralData.RawMassAxis);

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
        /// <exception cref="WorkspaceException">Thrown when some of the required data have not been loaded in the <see cref="Workspace"/> or their lengths are different.</exception>
        public void ResolutionFit(Interpolation.Type t, int order, double[] comList = null, double[] resolutionList = null)
        {
            // there are three use cases
            if (!IFDLoaded && comList == null && resolutionList == null) // this covers the use cases for IFJFile and GUI, when the COMList and ResolutionList have to be built.
            {
                if (Calibration.NameList.Count == 0) throw new WorkspaceException("Calibration namelist empty.");

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
            if (Calibration.COMList == null || Calibration.ResolutionList == null) throw new WorkspaceException("Resolution calibration points not specified.");
            if (Calibration.COMList.Length != Calibration.ResolutionList.Length) throw new WorkspaceException("Supplied resolution calibration point arrays have different lengths.");

            switch (t)
            {
                case Interpolation.Type.Polynomial:
                    Calibration.ResolutionInterp = new PolyInterpolation(Calibration.COMList, Calibration.ResolutionList, order);
                    break;
                case Interpolation.Type.SplineNatural:
                    throw new NotImplementedException("Natural spline interpolation has not yet been implemented.");
                //break;
                case Interpolation.Type.SplineNotAKnot:
                    Calibration.ResolutionInterp = new PPInterpolation(Calibration.COMList, Calibration.ResolutionList, PPInterpolation.PPType.SplineNotAKnot);
                    break;
                case Interpolation.Type.PCHIP:
                    Calibration.ResolutionInterp = new PPInterpolation(Calibration.COMList, Calibration.ResolutionList, PPInterpolation.PPType.PCHIP);
                    break;
                default:
                    throw new Interpolation.InterpolationException("Unknown interpolation type.");
            }
        }

        /// <summary>
        /// Builds the design matrix from currently supplied data and stores the result in the <see cref="Workspace.DesignMatrix"/> property.
        /// </summary>
        /// <remarks>The parameters <paramref name="searchRange"/> and <paramref name="fwhmRange"/> are being multiplied, so they must be both nonzero.</remarks>
        /// <param name="searchRange">Optional parameter, that sets the range around the centre of mass of a cluster, in which the design matrix values for that cluster will be generated.</param>
        /// <param name="fwhmRange">Optional parameter, that sets the range as a function of FWHM around the centre of mass of a cluster, in which the design matrix values for that cluster will be generated.</param>
        public void BuildDesignMatrix(double searchRange = 1, double fwhmRange = 0.5)
        {
            // TODO: merge this with the fit abundances maybe? If we do the design matrix updating, we will still need this function.
            DesignMatrix = new DesignMtrx(SpectralData, Clusters, Calibration)
            {
                SearchRange = searchRange,
                FwhmRange = fwhmRange
            };

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

            //Abundances = new OrderedDictionary(Clusters.Count);

            for (int i = 0; i < Clusters.Count; i++)
            {
                //Abundances.Add((Clusters[i] as IFData.Cluster).Name, lss.Solution[i]);
                (Clusters[i] as IFData.Cluster).Abundance = lss.Solution[i];
                //TODO: abundance error
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
