using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;

using MathNet.Numerics.LinearAlgebra;

using IsotopeFit;
using IsotopeFit.Numerics;

namespace IsotopeFitter
{
    class Program
    {
        //TODO: vymysliet, napisat, otestovat sposob rucneho zapisovania do datovych struktur

        static void Main(string[] args)
        {
            var bodka = new System.Globalization.NumberFormatInfo { NumberDecimalSeparator = "." };
            MathNet.Numerics.Control.UseNativeMKL();

            Stopwatch time = new Stopwatch();            

            Workspace W = new Workspace("testfiles\\finaltestfile.ifd");
            time.Start();

            W.CorrectBaseline();

            //TODO: the first fit needs to be spline (matlab spline(x,y,xx)) for that particular test file, not PCHIP. the second one is hardcoded pchip in matlab as well
            W.CorrectMassOffset(Interpolation.Type.SplineNotAKnot, 0);

            //TODO: resolutionfit

            Debug.Assert(W.SpectralData.MassAxis.Length == W.SpectralData.SignalAxis.Length);

            //using (FileStream f = File.Open("pureMassSig.txt", FileMode.Create))
            //{
            //    using (StreamWriter sw = new StreamWriter(f))
            //    {
            //        for (int i = 0; i < W.SpectralData.MassOffsetCorrAxis.Count; i++)
            //        {
            //            sw.WriteLine(W.SpectralData.MassOffsetCorrAxis[i].ToString(bodka) + " " + W.SpectralData.PureSignalAxis[i].ToString(bodka));
            //        }
            //    }
            //}

            W.BuildDesignMatrix();

            //using (FileStream f = File.Open("designMatrix.txt", FileMode.Create))
            //{
            //    using (StreamWriter sw = new StreamWriter(f))
            //    {
            //        int j = 0;

            //        for (int i = 0; i < W.designMatrix.Storage.Values.Length; i++)
            //        {
            //            if (i == W.designMatrix.Storage.ColumnPointers[j]) j++;
            //            sw.WriteLine(W.designMatrix.Storage.RowIndices[i].ToString(bodka) + " " + j.ToString(bodka) + " " + W.designMatrix.Storage.Values[i].ToString());
            //        }
            //    }
            //}

            W.FitAbundances();

            time.Stop();

            using (FileStream f = File.Open("abd.txt", FileMode.Create))
            {
                using (StreamWriter sw = new StreamWriter(f))
                {
                    for (int i = 0; i < W.Cluster.Count; i++)
                    {
                        sw.WriteLine(W.Cluster[i].CentreOfMass.ToString(bodka) + " " + W.Abundances[i].ToString(bodka));
                    }                    
                } 
            }

            //Vector<double> calcSpectrum = Vector<double>.Build.Dense((int)W.EndIndex);
            double[] calcSpectrum = new double[(int)W.EndIndex];
            //Vector<double> abdVec = Vector<double>.Build.Dense(W.Abundances);

            //now we have to remove the last vector from the design matrix

            double[] values = new double[W.DesignMatrix.Storage.ColumnPointers[W.DesignMatrix.Storage.ColumnPointers.Length - 2]];
            int[] rowIndices = new int[W.DesignMatrix.Storage.ColumnPointers[W.DesignMatrix.Storage.ColumnPointers.Length - 2]];
            int[] colPointers = new int[W.DesignMatrix.Storage.ColumnPointers.Length - 1];

            Array.Copy(W.DesignMatrix.Storage.Values, values, values.Length);
            Array.Copy(W.DesignMatrix.Storage.RowIndices, rowIndices, rowIndices.Length);
            Array.Copy(W.DesignMatrix.Storage.ColumnPointers, colPointers, colPointers.Length);

            CSparse.Double.SparseMatrix qwerty = new CSparse.Double.SparseMatrix(W.DesignMatrix.Storage.RowCount, W.DesignMatrix.Storage.ColumnCount - 1)
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

            

            using (FileStream f = File.Open("calcSpectrum.txt", FileMode.Create))
            {
                using (StreamWriter sw = new StreamWriter(f))
                {
                    for (int i = 0; i < W.SpectralData.RawLength; i++)
                    {
                        sw.WriteLine(W.SpectralData.MassAxis[i].ToString(bodka) + " " + calcSpectrum[i].ToString(bodka) + " " + W.SpectralData.SignalAxis[i].ToString(bodka));
                    }
                }
            }

            Console.WriteLine("done " + time.Elapsed);
            Console.ReadKey();
        }

    }
}
