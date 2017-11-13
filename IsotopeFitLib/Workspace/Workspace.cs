using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

namespace IsotopeFit
{
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

        #endregion

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

        //public void Dummy()
        //{
        //    DesignMtrx dm = new DesignMtrx(RawData.Length, Molecules);
        //    dm.Build();
        //}

        #endregion
    }
}
