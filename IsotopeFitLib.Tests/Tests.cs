using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Reflection;

using MathNet.Numerics.LinearAlgebra;

using IsotopeFit.Numerics;

namespace IsotopeFit.Tests
{
    [TestFixture]
    public partial class Tests
    {
        [Test, Category("Data manipulation")]
        public void BaselineSubtractTest()
        {
            Workspace Wrk = new Workspace(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\testfile.ifd");

            Wrk.CorrectBaseline();

            // solution check
            string[] bgCorrFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\bgcorrected.txt");

            List<double> bgCorr = new List<double>();

            foreach (string line in bgCorrFile)
            {
                if (!line.Contains("#") && line != "")
                {
                    bgCorr.Add(Convert.ToDouble(line.Trim(), new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                }
            }

            Assert.AreEqual(bgCorr.Count, Wrk.SpectralData.PureSignalAxis.Count);

            for (int i = 0; i < bgCorr.Count; i++)
            {
                Assert.AreEqual(bgCorr[i], Wrk.SpectralData.PureSignalAxis[i], 1e-9);
            }

            Assert.Pass("Baseline subtraction test passed");
        }

        [Test, Category("Data manipulation")]
        public void MassOffsetCorrectionTest()
        {
            Workspace Wrk = new Workspace(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\testfile.ifd");
            Wrk.CorrectBaseline();

            Wrk.CorrectMassOffset();

            // solution check
            string[] mOffFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\massAxisCorr.txt");

            List<double> mOff = new List<double>();

            foreach (string line in mOffFile)
            {
                if (!line.Contains("#") && line != "")
                {
                    mOff.Add(Convert.ToDouble(line.Trim(), new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                }
            }

            Assert.AreEqual(mOff.Count, Wrk.SpectralData.MassOffsetCorrAxis.Count);

            for (int i = 0; i < mOff.Count; i++)
            {
                Assert.AreEqual(mOff[i], Wrk.SpectralData.MassOffsetCorrAxis[i], 1e-9);
            }

            Assert.Pass("Mass offset correction test passed");
        }

        [Test, Category("Numerical algorithms")]
        public void NNLSTest()
        {
            //load test data from files
            string[] Afile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\nnls_C.txt");   //TODO: this can fail on Linux because of the backslashes
            string[] bfile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\nnls_d.txt");

            List<Vector<double>> Arows = new List<Vector<double>>();
            List<double> blist = new List<double>();

            foreach (string line in Afile)
            {
                List<double> values = new List<double>();

                if (!line.Contains("#") && line != "")
                {
                    string[] valuesstr = line.Trim().Split(new char[] { ' ' });

                    foreach (string str in valuesstr)
                    {
                        values.Add(Convert.ToDouble(str, new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                    }

                    Arows.Add(Vector<double>.Build.DenseOfEnumerable(values));
                }
            }

            foreach (string line in bfile)
            {
                List<double> values = new List<double>();

                if (!line.Contains("#") && line != "")
                {
                    blist.Add(Convert.ToDouble(line, new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                }
            }

            Matrix<double> C = Matrix<double>.Build.DenseOfRowVectors(Arows);
            Vector<double> d = Vector<double>.Build.DenseOfEnumerable(blist);

            Workspace wrk = new Workspace();
            LeastSquaresSystem lss = new LeastSquaresSystem(C, d);
            lss.Solve();
            Vector<double> solution = lss.Solution;

            // solution check
            string[] xfile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\nnls_x_solution.txt");

            List<double> correctX = new List<double>();

            foreach (string line in xfile)
            {
                correctX.Add(Convert.ToDouble(line));   //TODO: this might fail if decimal symbol is wrong
            }

            for (int i = 0; i < correctX.Count; i++)
            {
                //Assert.Less(1e-9, Math.Abs(solution[i] - correctX[i]),  "Solution is wrong.");
                Assert.AreEqual(correctX[i], solution[i], 1e-9);
            }

            Assert.Pass("NNLS test passed.");
        }

        [Test]
        public void IFDLoadTest()
        {
            Workspace Wrk = new Workspace(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\testfile.ifd");

            // raw data check
            Assert.AreEqual(4.64416598796176d, Wrk.SpectralData.RawMassAxis[9999], 1e-12);
            Assert.AreEqual(173.00000000000000d, Wrk.SpectralData.RawSignalAxis[9999], 1e-12);

            // start and end index check
            Assert.AreEqual(1, Wrk.StartIndex);
            Assert.AreEqual(416375, Wrk.EndIndex);

            // molecules check
            IFData.Molecule M = Wrk.Molecules[121];
            Assert.AreEqual(2.45148615736600e+2d, M.PeakData.Mass[1], 1e-12);
            Assert.AreEqual(8.80748376064873e-2d, M.PeakData.Abundance[1], 1e-12);
            Assert.AreEqual("[TMAl][DMAl - 1H]3[H]4", M.Name);
            Assert.AreEqual(244.145260901600d, M.MinMass, 1e-12);
            Assert.AreEqual(246.151970571600d, M.MaxMass, 1e-12);
            Assert.AreEqual(67539, M.MinIndex);
            Assert.AreEqual(67813, M.MaxIndex);
            Assert.AreEqual(2683.29332483441d, M.Area, 1e-10);
            Assert.AreEqual(1111.36756079838d, M.AreaError, 1e-10);
            Assert.AreEqual(122, M.RootIndex);

            // calibration check
            Assert.AreEqual(168.109336672200d, Wrk.Calibration.COMList[24], 1e-12);
            Assert.AreEqual(-0.170476870963361d, Wrk.Calibration.MassOffsetList[38], 1e-12);
            Assert.AreEqual(2711.93900551692d, Wrk.Calibration.ResolutionList[12], 1e-12);
            Assert.AreEqual("Spline", Wrk.Calibration.MassOffsetMethod);
            Assert.AreEqual("Polynomial", Wrk.Calibration.ResolutionMethod);
            Assert.AreEqual(0, Wrk.Calibration.MassOffsetParam);
            Assert.AreEqual(2, Wrk.Calibration.ResolutionParam);
            Assert.AreEqual("[DMAl][N2][TMAl][H]", Wrk.Calibration.Namelist[46]);

            //shape check
            Assert.AreEqual("pp", Wrk.Calibration.Shape.Form);
            Assert.AreEqual(0.0867891545686533d, Wrk.Calibration.Shape.Breaks[8], 1e-12);
            Assert.AreEqual(-0.07977458312173122d, Wrk.Calibration.Shape.Coefs[2, 1], 1e-12);
            Assert.AreEqual(15, Wrk.Calibration.Shape.Pieces);
            Assert.AreEqual(4, Wrk.Calibration.Shape.Order);
            Assert.AreEqual(1, Wrk.Calibration.Shape.Dim);

            // bgcorrection check
            Assert.AreEqual(Double.NegativeInfinity, Wrk.BaselineCorrData.StartMass, 1e-12);
            Assert.AreEqual(Double.PositiveInfinity, Wrk.BaselineCorrData.EndMass, 1e-12);
            Assert.AreEqual(50, Wrk.BaselineCorrData.NDiv);
            Assert.AreEqual(40d, Wrk.BaselineCorrData.Percent);
            Assert.AreEqual(2450.15730478338d, Wrk.BaselineCorrData.XAxis[25], 1e-10);
            Assert.AreEqual(190.835184629240d, Wrk.BaselineCorrData.YAxis[42], 1e-10);

            Assert.Pass("IFDFile read test passed.");
        }
    }
}
