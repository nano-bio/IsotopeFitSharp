using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Reflection;

using CSparse.Double;

namespace IsotopeFitLib.Tests
{
    [TestFixture]
    public class MatrixInversionTests
    {
        [Test, Category("Numerical algorithms")]
        public void MatrixInversion1()
        {
            string[] valFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\MatrixInversion\\mtival1.txt");   //TODO: this can fail on Linux because of the backslashes
            string[] ridxFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\MatrixInversion\\mtiridx1.txt");
            string[] cptFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\MatrixInversion\\mticp1.txt");

            List<double> values = new List<double>();
            List<int> rowInd = new List<int>();
            List<int> colPt = new List<int>();

            for (int i = 0; i < valFile.Length; i++)
            {
                values.Add(Convert.ToDouble(valFile[i], new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                rowInd.Add(Convert.ToInt32(ridxFile[i]));                
            }

            for (int i = 0; i < cptFile.Length; i++)
            {
                colPt.Add(Convert.ToInt32(cptFile[i]));
            }

            SparseMatrix A = new SparseMatrix(colPt.Count - 1, colPt.Count - 1)
            {
                Values = values.ToArray(),
                RowIndices = rowInd.ToArray(),
                ColumnPointers = colPt.ToArray()
            };

            // TODO: call the matrix inverse
            SparseMatrix In = IsotopeFit.MatrixInversion.Inverse(A);

            string[] ivalFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\MatrixInversion\\imtval1.txt");   //TODO: this can fail on Linux because of the backslashes
            string[] iridxFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\MatrixInversion\\imtridx1.txt");
            string[] icptFile = File.ReadAllLines(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\MatrixInversion\\imtcpt1.txt");

            List<double> ivalues = new List<double>();
            List<int> irowInd = new List<int>();
            List<int> icolPt = new List<int>();

            for (int i = 0; i < ivalFile.Length; i++)
            {
                ivalues.Add(Convert.ToDouble(ivalFile[i], new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." }));
                irowInd.Add(Convert.ToInt32(iridxFile[i]));
            }

            for (int i = 0; i < icptFile.Length; i++)
            {
                icolPt.Add(Convert.ToInt32(icptFile[i]));
            }

            // Assertions block
            for (int i = 0; i < ivalues.Count; i++)
            {
                Assert.AreEqual(ivalues[i], In.Values[i], 1e-6);
                Assert.AreEqual(irowInd[i], In.RowIndices[i]);                
            }

            for (int i = 0; i < icolPt.Count; i++)
            {
                Assert.AreEqual(icolPt[i], In.ColumnPointers[i]);
            }

            Assert.Pass("Matrix inversion test 1 passed.");
        }
    }
}
