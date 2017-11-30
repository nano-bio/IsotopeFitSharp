using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using System.Reflection;

using CSparse.Double;

using IsotopeFit;
using IsotopeFit.Numerics;

namespace IsotopeFitter.Tests
{
    [TestFixture]
    public class Tests
    {
        System.Globalization.NumberFormatInfo dot = new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." };

        [Test, Category("IsotopeFitter")]
        public void WorkflowTests()
        {
            MathNet.Numerics.Control.UseNativeMKL();

            Workspace W = new Workspace(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\finaltestfile.ifd");

            Stopwatch timeTotal = new Stopwatch();
            timeTotal.Start();

            BaselineSubtest(ref W);

            MassOffsetSubtest(ref W);

            //TODO: resolution fit subtest
            ResolutionFitSubtest(ref W);

            DesignMatrixBuildSubtest(ref W);            

            ExtractAbundancesSubtest(ref W);
            
            timeTotal.Stop();
            Assert.Pass("Workflow test passed. Elapsed time: {0}", timeTotal.Elapsed);
        }

        private void BaselineSubtest(ref Workspace w)
        {
            w.CorrectBaseline();

            // compare calculated pure signal with matlab results
            string[] bgCorrFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\1outSubtractBg.txt");

            List<double> bgCorr = new List<double>();

            foreach (string line in bgCorrFile)
            {
                if (!line.Contains("#") && line != "")
                {
                    bgCorr.Add(Convert.ToDouble(line.Trim(), dot));
                }
            }

            Assert.AreEqual(bgCorr.Count, w.SpectralData.SignalAxis.Length);

            for (int i = 0; i < bgCorr.Count; i++)
            {
                Assert.AreEqual(bgCorr[i],w.SpectralData.SignalAxis[i], 1e-9, "baseline check failed at index {0}", i);
            }
        }


        private void MassOffsetSubtest(ref Workspace w)
        {
            //TODO: the first fit needs to be spline (matlab spline(x,y,xx)) for that particular test file, not PCHIP. the second one is hardcoded pchip in matlab as well
            w.CorrectMassOffset(Interpolation.Type.SplineNotAKnot, 0);

            // compare calculated mass axis with matlab results
            string[] massOffFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\2outSubtractMassOffset.txt");

            List<double> mOff = new List<double>();

            foreach (string line in massOffFile)
            {
                if (!line.Contains("#") && line != "")
                {
                    mOff.Add(Convert.ToDouble(line.Trim(), dot));
                }
            }

            Assert.AreEqual(mOff.Count, w.SpectralData.MassAxis.Length);

            for (int i = 0; i < mOff.Count; i++)
            {
                Assert.AreEqual(mOff[i], w.SpectralData.MassAxis[i], 1e-9, "mass offset check failed at index {0}", i);
            }
        }


        private void ResolutionFitSubtest(ref Workspace w)
        {
            w.ResolutionFit(Interpolation.Type.Polynomial, 2);

            // compare calculated resolution coefficients with matlab results
            string[] resCoefFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\3resolutionCoefs.txt");

            List<double> resCoef = new List<double>();

            foreach (string line in resCoefFile)
            {
                if (!line.Contains("#") && line != "")
                {
                    resCoef.Add(Convert.ToDouble(line.Trim(), dot));
                }
            }

            resCoef.Reverse();

            Assert.AreEqual(resCoef.Count, (w.ResolutionInterpolation as PolyInterpolation).Coefs.Length);

            for (int i = 0; i < resCoef.Count; i++)
            {
                Assert.AreEqual(resCoef[i], (w.ResolutionInterpolation as PolyInterpolation).Coefs[i], 1e-9, "mass offset check failed at index {0}", i);
            }
        }

        private void DesignMatrixBuildSubtest(ref Workspace w)
        {
            w.BuildDesignMatrix();

            // compare with matlab calculated abundances
            string[] matrixFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\5designMatrixWithMask.txt");

            double[] values = new double[matrixFile.Length];
            int[] rowIndices = new int[matrixFile.Length];
            int[] colIndices = new int[matrixFile.Length];

            for (int i = 0; i < matrixFile.Length; i++)
            {
                string[] hue = matrixFile[i].Trim().Split(' ');

                values[i] = Convert.ToDouble(hue[2], dot);
                rowIndices[i] = (int)Convert.ToDouble(hue[0], dot) - 1;
                colIndices[i] = (int)Convert.ToDouble(hue[1], dot) - 1;
            }

            int cols = colIndices.Distinct().Count();
            int[] colPointers = new int[cols + 1];
            colPointers[0] = 0;

            var colCounts = colIndices.GroupBy(x => x).ToArray();
            int cumsum = 0;

            foreach (var item in colCounts)
            {
                cumsum += item.Count();
                colPointers[item.Key + 1] = cumsum;
            }

            int bue = w.DesignMatrix.Storage.RowIndices.Max();

            // comparisons of the matrix internal data fields
            Assert.AreEqual(values.Length, w.DesignMatrix.Storage.Values.Length - (w.DesignMatrix.Storage.ColumnPointers.Last() - w.DesignMatrix.Storage.ColumnPointers[w.DesignMatrix.Storage.ColumnPointers.Length - 2]), "different values array length in the design matrix"); //TODO: this will be a problem with non-zero observation vector
            Assert.AreEqual(rowIndices.Length, w.DesignMatrix.Storage.RowIndices.Length - (w.DesignMatrix.Storage.ColumnPointers.Last() - w.DesignMatrix.Storage.ColumnPointers[w.DesignMatrix.Storage.ColumnPointers.Length - 2]), "different row indices array length in the design matrix");
            Assert.AreEqual(colPointers.Length, w.DesignMatrix.Storage.ColumnPointers.Length - 1, "different column pointers array length in the design matrix");   // -1 because the calculated design matrix contains an extra column for the observation vector

            for (int i = 0; i < rowIndices.Length; i++)
            {
                Assert.AreEqual(rowIndices[i], w.DesignMatrix.Storage.RowIndices[i], "design matrix row index comparison fail at {0}", i);
                Assert.AreEqual(values[i], w.DesignMatrix.Storage.Values[i], 1e-9, "design matrix values comparison fail at {0}", i);
            }

            for (int i = 0; i < colPointers.Length; i++)
            {
                Assert.AreEqual(colPointers[i], w.DesignMatrix.Storage.ColumnPointers[i], "design matrix column pointers comparison fail at {0}", i);
            }

            // check the observation vector
            string[] obsVecFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\5signalWithMask.txt");

            List<double> obsVecCorrect = new List<double>();

            foreach (string line in obsVecFile)
            {
                if (!line.Contains("#") && line != "")
                {
                    obsVecCorrect.Add(Convert.ToDouble(line.Trim(), dot));
                }
            }

            int lastColumnCount = w.DesignMatrix.Storage.ColumnPointers.Last() - w.DesignMatrix.Storage.ColumnPointers[w.DesignMatrix.Storage.ColumnPointers.Length - 2];
            double[] obsVec = new double[lastColumnCount];

            Array.Copy(w.DesignMatrix.Storage.Values, w.DesignMatrix.Storage.Values.Length - lastColumnCount, obsVec, 0, lastColumnCount);

            for (int i = 0; i < obsVecCorrect.Count; i++)
            {
                Assert.AreEqual(obsVecCorrect[i], obsVec[i], 1e-6, "observation vector comparison fail at {0}", i);
            }
        }

        private void ExtractAbundancesSubtest(ref Workspace w)
        {
            w.FitAbundances();

            // compare with matlab calculated abundances
            string[] abdFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\6abundancesFromLsqnonneg.txt");

            List<double> abd = new List<double>();

            foreach (string line in abdFile)
            {
                if (!line.Contains("#") && line != "")
                {
                    abd.Add(Convert.ToDouble(line.Trim(), dot));
                }
            }

            Assert.AreEqual(abd.Count, w.Abundances.Length);

            for (int i = 0; i < abd.Count; i++)
            {
                Assert.AreEqual(abd[i], w.Abundances[i], 1e-6, "abundances check failed at index {0}", i);
            }
        }
    }
}
