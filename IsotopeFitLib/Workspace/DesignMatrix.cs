using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

using MathNet.Numerics.LinearAlgebra;

namespace IsotopeFit
{
    public partial class Workspace
    {
        internal class DesignMtrx
        {
            //int peakLowerLimitIndex, peakUpperLimitIndex, fitmaskLowerLimitIndex, fitmaskUpperLimitIndex, breakIndex;
            //double mass, area, resolution, fwhm, idealSignalValue;

            //bool[] fitMask;

            /// <summary>
            /// Create new empty design matrix.
            /// </summary>
            internal DesignMtrx()
            {
                buildColumnDelegate = BuildColumn;
            }

            internal Matrix<double> Storage { get; set; }
            internal Matrix<double> R { get; set; }

            //private delegate Vector<double> BuildColumnDelegate(int index);
            private Action<int> buildColumnDelegate;

            internal void Build()
            {
                //TODO: we will probably want to use the parallel for overload with the init-body-final scheme
                Parallel.For(0, 100, buildColumnDelegate);
                //Parallel.For(0, 100,)
            }

            internal void CalculateQR()
            {
                throw new NotImplementedException();
            }

            private void BuildColumn(int moleculeIndex)
            {
                Console.WriteLine(Thread.CurrentThread.ManagedThreadId);
                
            }
        }
    }
}
