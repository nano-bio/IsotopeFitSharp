using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Reflection;

using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Storage;

using IsotopeFit.Numerics;

namespace IsotopeFitLib.Tests
{
    [TestFixture]
    public partial class Tests
    {
        [Test]
        public void LeaSqrQRTest()
        {
            //load test data from files
            string[] Afile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\designMatrix.txt");   //TODO: this can fail on Linux because of the backslashes
            string[] bfile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\spectrum.txt");

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

            // original matrix
            SparseMatrix C = (SparseMatrix)SparseMatrix.Build.SparseOfRowVectors(Arows);
            Vector<double> d = Vector<double>.Build.DenseOfEnumerable(blist);

            LeastSquaresSystem lss = new LeastSquaresSystem(C, d);

            lss = Algorithm.LeaSqrSparseQR(lss);            

            // solve the system
            lss.Solve();

            // solution check
            string[] xfile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\abundances.txt");

            List<double> correctX = new List<double>();
            
            foreach (string line in xfile)
            {
                if (!line.Contains("#") && line != "")
                {
                    correctX.Add(Convert.ToDouble(line.Trim(), new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                }
            }

            for (int i = 0; i < correctX.Count; i++)
            {
                //Assert.Less(1e-9, Math.Abs(solution[i] - correctX[i]),  "Solution is wrong.");
                Assert.AreEqual(correctX[i], lss.Solution[i], 1e-9);
            }

            Assert.Pass("NNLS with sparse QR factorization test passed.");
        }
    }
}
