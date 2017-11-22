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
        private DesignMtrx DesignMatrix;

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

        public IFData.Spectrum RawData { get; set; }
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

            RawData = IFDFile.ReadRawData(rootElement);
            StartIndex = IFDFile.ReadStartIndex(rootElement);
            EndIndex = IFDFile.ReadEndIndex(rootElement);
            Molecules = IFDFile.ReadMolecules(rootElement);
            Calibration = IFDFile.ReadCalibration(rootElement);
            BaselineCorr = IFDFile.ReadBackgroundCorr(rootElement);
        }


        public void CorrectBaseline()
        {
            //TODO: check is all necessary data are already loaded in the workspace. throw exception if not

            //TODO: calculate the background corrected spectrum

            //TODO: store it in a corresponding workspace field              
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
        }

        public void BuildDesignMatrix()
        {
            DesignMtrx dm = new DesignMtrx(RawData, Molecules, Calibration);
            dm.Build();
        }

        #endregion
    }
}
