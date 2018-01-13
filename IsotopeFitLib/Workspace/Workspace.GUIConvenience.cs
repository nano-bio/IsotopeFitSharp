using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IsotopeFit
{
    public partial class Workspace
    {
        /// <summary>
        /// Convenience method for the GUI that calculates estimate of a spectrum fit from new mass offset and resolution data and old abundance values.
        /// </summary>
        /// <param name="massOffsetInterpType">Type of the mass offset interpolation.</param>
        /// <param name="resInterpType">Type of the resolution interpolation.</param>
        /// <param name="resInterpOrder">Order of the resolution interpolation, if polynomial is used. Otherwise ignored.</param>
        public void GUICalibEstimateSpectrumFit(Interpolation.Type massOffsetInterpType, Interpolation.Type resInterpType, int massOffsetInterpOrder = -1, int resInterpOrder = -1, bool massAxisAutoCrop = false)
        {
            CorrectMassOffset(massOffsetInterpType, massOffsetInterpOrder, massAxisAutoCrop);
            ResolutionFit(resInterpType, resInterpOrder);
            BuildDesignMatrix();    //TODO: the fwhmRange and searchRange should also be settable, either here, or in some more central way
            CalculateSpectrum();
        }

        
    }
}
