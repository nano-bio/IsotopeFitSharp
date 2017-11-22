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
        [Test]
        public void ResolutionPolynomialFitTest()
        {
            Workspace Wrk = new Workspace(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\testfile.ifd");

            PolyInterpolation PolyRC = new PolyInterpolation(Wrk.Calibration.COMList.ToArray(), Wrk.Calibration.ResolutionList.ToArray(), 3);

            double[] Solution = PolyRC.Coefs;
            Solution = Solution.Reverse().ToArray();

            string[] PolyCoefs = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\PolynomialResolutionCoefs.txt");

            List<double> CorrectedCoefs = new List<double>();

            foreach (string line in PolyCoefs)
            {
                if (!line.Contains("#") && line != "")
                {
                    string[] valuesstr = line.Trim().Split(new char[] { ' ' });

                    foreach (string str in valuesstr)
                    {
                        CorrectedCoefs.Add(Convert.ToDouble(str, new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                    }
                }
            }

            for (int i = 0; i < CorrectedCoefs.ToArray().Length; i++)
            {
                Assert.AreEqual(CorrectedCoefs.ToArray()[i], Solution[i], 1e-9);
            }

            Assert.Pass("Polynomial resolution fit test passed.");
        }

        [Test]
        public void ResolutionPCHIPFitTest()
        {
            Workspace Wrk = new Workspace(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\testfile.ifd");

            PPInterpolation PPRC = new PPInterpolation(Wrk.Calibration.COMList.ToArray(), Wrk.Calibration.ResolutionList.ToArray(), PPInterpolation.PPType.PCHIP);

            double[][] Solution = PPRC.Coefs;
            string[] PCHIPCoefs = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\PCHIPResolutionCoefs.txt");

            List<Vector<double>> rows = new List<Vector<double>>();

            foreach (string line in PCHIPCoefs)
            {
                List<double> values = new List<double>();

                if (!line.Contains("#") && line != "")
                {
                    string[] valuesstr = line.Trim().Split(new char[] { ' ' });

                    foreach (string str in valuesstr)
                    {
                        values.Add(Convert.ToDouble(str, new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                    }

                    rows.Add(Vector<double>.Build.DenseOfEnumerable(values));
                }
            }

            double[][] CorrectedCoefs = Matrix<double>.Build.DenseOfRowVectors(rows).ToRowArrays();
            
            Assert.Pass("PCHIP resolution fit test passed.");
        }

        [Test]
        public void EvaluationPolynomialTest()
        {
            Workspace Wrk = new Workspace(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\testfile.ifd");

            PolyInterpolation PolyRC = new PolyInterpolation(Wrk.Calibration.COMList.ToArray(), Wrk.Calibration.ResolutionList.ToArray(), 3);

            double[] Solution = PolyRC.Evaluate(Wrk.Calibration.COMList.ToArray());

            string[] PolyEval = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\PolynomialEvaluation.txt");

            List<double> CorrectedEval = new List<double>();

            foreach (string line in PolyEval)
            {
                if (!line.Contains("#") && line != "")
                {
                    string[] valuesstr = line.Trim().Split(new char[] { ' ' });

                    foreach (string str in valuesstr)
                    {
                        CorrectedEval.Add(Convert.ToDouble(str, new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                    }
                }
            }

            for (int i = 0; i < CorrectedEval.ToArray().Length; i++)
            {
                Assert.AreEqual(CorrectedEval.ToArray()[i], Solution[i], 1e-9);
            }

            Assert.Pass("Polynomial Evaluation test passed.");
        }

        [Test]
        public void EvaluationPCHIPTest()
        {
            Workspace Wrk = new Workspace(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\testfile.ifd");

            PPInterpolation PPRC = new PPInterpolation(Wrk.Calibration.COMList.ToArray(), Wrk.Calibration.ResolutionList.ToArray(), PPInterpolation.PPType.PCHIP);

            string[] PCHIPCoefs = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\PCHIPEvaluation.txt");

            List<Vector<double>> rows = new List<Vector<double>>();

            foreach (string line in PCHIPCoefs)
            {
                List<double> values = new List<double>();

                if (!line.Contains("#") && line != "")
                {
                    string[] valuesstr = line.Trim().Split(new char[] { ' ' });

                    foreach (string str in valuesstr)
                    {
                        values.Add(Convert.ToDouble(str, new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                    }

                    rows.Add(Vector<double>.Build.DenseOfEnumerable(values));
                }
            }
            
            Matrix<double> M = Matrix<double>.Build.DenseOfRowVectors(rows);
            
            double[] Solution = PPRC.Evaluate(M.Column(0).ToArray());

            for (int i = 0; i < Solution.GetLength(0); i++)
            {
                Assert.AreEqual(Convert.ToDouble(M.Column(1).ToArray()[i], new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }), Solution[i], 1e-9);
            }

            Assert.Pass("PCHIP Evaluation test passed.");
        }
    }
}
