using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;

using MathNet.Numerics.LinearAlgebra;

using IsotopeFit;

namespace IsotopeFitter
{
    class Program
    {
        //TODO: vymysliet, napisat, otestovat sposob rucneho zapisovania do datovych struktur

        static void Main(string[] args)
        {
            MathNet.Numerics.Control.UseNativeMKL();

            Stopwatch time = new Stopwatch();            

            Workspace W = new Workspace("testfiles\\finaltestfile.ifd");
            time.Start();

            W.CorrectBaseline();

            //TODO: the first fit needs to be spline (matlab spline(x,y,xx)) for that particular test file, not PCHIP. the second one is hardcoded pchip in matlab as well
            W.CorrectMassOffset();

            W.BuildDesignMatrix();

            W.ExtractAbundances();

            time.Stop();

            using (FileStream f = File.OpenWrite("abd.txt"))
            {
                using (StreamWriter sw = new StreamWriter(f))
                {
                    for (int i = 0; i < W.Molecules.Count; i++)
                    {
                        sw.WriteLine(W.Molecules[i].CentreOfMass.ToString() + " " + W.Abundances[i].ToString());
                    }                    
                } 
            }

            //Vector<double> calcSpectrum = Vector<double>.Build.Dense((int)W.EndIndex);
            double[] calcSpectrum = new double[(int)W.EndIndex];
            //Vector<double> abdVec = Vector<double>.Build.Dense(W.Abundances);

            //now we have to remove the last vector from the design matrix

            double[] values = new double[W.designMatrix.Storage.ColumnPointers[W.designMatrix.Storage.ColumnPointers.Length - 2]];
            int[] rowIndices = new int[W.designMatrix.Storage.ColumnPointers[W.designMatrix.Storage.ColumnPointers.Length - 2]];
            int[] colPointers = new int[W.designMatrix.Storage.ColumnPointers.Length - 1];

            Array.Copy(W.designMatrix.Storage.Values, values, values.Length);
            Array.Copy(W.designMatrix.Storage.RowIndices, rowIndices, rowIndices.Length);
            Array.Copy(W.designMatrix.Storage.ColumnPointers, colPointers, colPointers.Length);

            CSparse.Double.SparseMatrix qwerty = new CSparse.Double.SparseMatrix(W.designMatrix.Storage.RowCount, W.designMatrix.Storage.ColumnCount - 1)
            {
                Values = values,
                RowIndices = rowIndices,
                ColumnPointers = colPointers
            };

            //Matrix<double> mehehe = W.designMatrix.Storage.SubMatrix(0, W.designMatrix.Storage.RowCount, 0, W.designMatrix.Storage.ColumnCount - 1);

            //if (W.Abundances.Sum() == 0)
            //{
            //    Console.WriteLine("abd vsetko nuly");
            //}

            qwerty.Multiply(W.Abundances, calcSpectrum);

            //Vector<double> calcSpectrumVect = Vector<double>.Build.SparseOfArray(calcSpectrum);

            //if (calcSpectrum.Sum() == 0)
            //{
            //    Console.WriteLine("calc spec vsetko nuly");
            //}

            //calcSpectrum = qwerty.Multiply(abdVec);

            var bodka = new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." };

            using (FileStream f = File.OpenWrite("calcSpectrum.txt"))
            {
                using (StreamWriter sw = new StreamWriter(f))
                {
                    for (int i = 0; i < W.SpectralData.RawLength; i++)
                    {
                        sw.WriteLine(W.SpectralData.MassOffsetCorrAxis[i].ToString(bodka) + " " + calcSpectrum[i].ToString(bodka) + " " + W.SpectralData.PureSignalAxis[i].ToString(bodka));
                    }
                }
            }

            Console.WriteLine("done " + time.Elapsed);
            Console.ReadKey();
        }

    }
}
