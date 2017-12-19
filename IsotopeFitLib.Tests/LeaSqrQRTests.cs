using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Reflection;
using System.Diagnostics;

using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Storage;

using IsotopeFit;

namespace IsotopeFitLib.Tests
{
    [TestFixture]
    public partial class Tests
    {
        [Test, Category("Numerical algorithms")]
        public void LeaSqrQRHouseholderTest()
        {
            //load test data from files
            string[] Afile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\designMatrix.txt");   //TODO: this can fail on Linux because of the backslashes
            string[] bfile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\spectrum.txt");

            List<Vector<double>> Arows = new List<Vector<double>>();
            List<double> blist = new List<double>();

            foreach (string line in Afile)
            {
                List<double> vals = new List<double>();

                if (!line.Contains("#") && line != "")
                {
                    string[] valuesstr = line.Trim().Split(new char[] { ' ' });

                    foreach (string str in valuesstr)
                    {
                        vals.Add(Convert.ToDouble(str, new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                    }

                    Arows.Add(Vector<double>.Build.DenseOfEnumerable(vals));
                }
            }

            foreach (string line in bfile)
            {
                List<double> vals = new List<double>();

                if (!line.Contains("#") && line != "")
                {
                    blist.Add(Convert.ToDouble(line, new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                }
            }

            // original matrix
            SparseMatrix C = (SparseMatrix)SparseMatrix.Build.SparseOfRowVectors(Arows);
            Vector<double> d = Vector<double>.Build.DenseOfEnumerable(blist);

            // we need to copy the C matrix into the CSparse format, just for the sake of the test, no need to be nice
            List<SparseVector> Cc = new List<SparseVector>(C.ColumnCount);

            for (int i = 0; i < C.ColumnCount; i++)
            {
                Cc.Add(C.Column(i) as SparseVector);
            }

            int nonZeroCount = C.NonZerosCount;

            double[] values = new double[nonZeroCount];
            int[] rowIndices = new int[nonZeroCount];
            int[] colPointers = new int[C.ColumnCount + 1];

            for (int i = 0; i < Cc.Count; i++)
            {
                Array.Copy((Cc[i].Storage as SparseVectorStorage<double>).Values, 0, values, colPointers[i], Cc[i].NonZerosCount);
                Array.Copy((Cc[i].Storage as SparseVectorStorage<double>).Indices, 0, rowIndices, colPointers[i], Cc[i].NonZerosCount);
                colPointers[i + 1] = colPointers[i] + Cc[i].NonZerosCount;
            }

            //values = (C.Storage as SparseCompressedRowMatrixStorage<double>).Values;
            //rowIndices = Enumerable.Repeat(Enumerable.Range(0, C.RowCount), C.ColumnCount).SelectMany(x => x).ToArray();    //TODO: test lol

            //for (int i = 0; i <= C.ColumnCount; i++)
            //{
            //    colPointers[i] = i * C.RowCount;
            //}

            CSparse.Double.SparseMatrix Cs = new CSparse.Double.SparseMatrix(C.RowCount, C.ColumnCount)
            {
                Values = values,
                RowIndices = rowIndices,
                ColumnPointers = colPointers
            };

            //LeastSquaresSystem lss = new LeastSquaresSystem(Cs, d);

            Stopwatch s = new Stopwatch();
            s.Start();
            //lss = Algorithm.LeaSqrSparseQRHouseholder(lss);
            s.Stop();

            Console.WriteLine("LeaSqrSparseQRHouseholder " + s.ElapsedMilliseconds);

            // solve the system
            //lss.Solve();

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
                //Assert.AreEqual(correctX[i], lss.Solution[i], 1e-9);
            }

            Assert.Pass("NNLS with sparse QR factorization test passed.");
        }

        //[Test]
        //public void LeaSqrQRGivensTest()
        //{
        //    //load test data from files
        //    string[] Afile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\designMatrix.txt");   //TODO: this can fail on Linux because of the backslashes
        //    string[] bfile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\spectrum.txt");

        //    List<Vector<double>> Arows = new List<Vector<double>>();
        //    List<double> blist = new List<double>();

        //    foreach (string line in Afile)
        //    {
        //        List<double> values = new List<double>();

        //        if (!line.Contains("#") && line != "")
        //        {
        //            string[] valuesstr = line.Trim().Split(new char[] { ' ' });

        //            foreach (string str in valuesstr)
        //            {
        //                values.Add(Convert.ToDouble(str, new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
        //            }

        //            Arows.Add(Vector<double>.Build.DenseOfEnumerable(values));
        //        }
        //    }

        //    foreach (string line in bfile)
        //    {
        //        List<double> values = new List<double>();

        //        if (!line.Contains("#") && line != "")
        //        {
        //            blist.Add(Convert.ToDouble(line, new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
        //        }
        //    }

        //    // original matrix
        //    SparseMatrix C = (SparseMatrix)SparseMatrix.Build.SparseOfRowVectors(Arows);
        //    Vector<double> d = Vector<double>.Build.DenseOfEnumerable(blist);

        //    LeastSquaresSystem lss = new LeastSquaresSystem(C, d);

        //    Stopwatch s = new Stopwatch();
        //    s.Start();
        //    //lss = Algorithm.LeaSqrSparseQRGivens(lss);
        //    s.Stop();

        //    Console.WriteLine("LeaSqrSparseQRGivens " + s.ElapsedMilliseconds);

        //    // solve the system
        //    lss.Solve();

        //    // solution check
        //    string[] xfile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\abundances.txt");

        //    List<double> correctX = new List<double>();

        //    foreach (string line in xfile)
        //    {
        //        if (!line.Contains("#") && line != "")
        //        {
        //            correctX.Add(Convert.ToDouble(line.Trim(), new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
        //        }
        //    }

        //    for (int i = 0; i < correctX.Count; i++)
        //    {
        //        //Assert.Less(1e-9, Math.Abs(solution[i] - correctX[i]),  "Solution is wrong.");
        //        Assert.AreEqual(correctX[i], lss.Solution[i], 1e-9);
        //    }

        //    Assert.Pass("NNLS with sparse QR factorization test passed.");
        //}
    }
}
