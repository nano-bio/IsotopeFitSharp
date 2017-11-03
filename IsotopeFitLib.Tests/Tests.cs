using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Reflection;

using MathNet.Numerics.LinearAlgebra;

namespace IsotopeFitLib.Tests
{
    [TestFixture]
    public class Tests
    {
        [Test]
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

            IsotopeFitWorkspace wrk = new IsotopeFitWorkspace();
            Vector<double> solution = wrk.NNLS(C, d);

            // solution check
            string[] xfile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\nnls_x_solution.txt");

            List<double> correctX = new List<double>();

            foreach (string line in xfile)
            {
                correctX.Add(Convert.ToDouble(line));   //TODO: this might fail if decimal symbol is wrong
            }

            for (int i = 0; i < correctX.Count; i++)
            {
                Assert.Less(Math.Abs(solution[i] - correctX[i]), 1e-9, "Solution is wrong.");
            }

            Assert.Pass("NNLS test passed.");
        }
    }
}
