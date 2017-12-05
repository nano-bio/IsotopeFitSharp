using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using System.Globalization;

using MathNet.Numerics.LinearAlgebra;

using IsotopeFit;
//using IsotopeFit.Numerics;

namespace IsotopeFitter
{
    class Program
    {
        NumberFormatInfo dot = new NumberFormatInfo { NumberDecimalSeparator = "." };

        static void Main(string[] args)
        {
            if (args.Length == 0)   //TODO: this might fail on linux
            {
                throw new Exception("No file or directory to process was specified.");
            }

            Stopwatch baselineTime = new Stopwatch();
            Stopwatch mofTime = new Stopwatch();
            Stopwatch resTime = new Stopwatch();
            Stopwatch dmTime = new Stopwatch();
            Stopwatch abdTime = new Stopwatch();

            string command = args[0];

            Workspace W = new Workspace();

            if (command.Contains(".ifd"))
            {
                W.LoadIFDFile(command);
            }
            else if (command.Contains(".ifj"))
            {
                throw new NotImplementedException();
            }
            else
            {
                // command is a directory
                throw new NotImplementedException();
            }


            
            W.CorrectBaseline();
            baselineTime.Stop();
            Console.WriteLine("baseline done in {0} miliseconds", baselineTime.ElapsedMilliseconds);

            mofTime.Start();
            W.CorrectMassOffset(Interpolation.Type.SplineNotAKnot);
            mofTime.Stop();
            Console.WriteLine("mass offset done in {0} miliseconds", mofTime.ElapsedMilliseconds);

            resTime.Start();
            W.ResolutionFit(Interpolation.Type.Polynomial, 2);
            resTime.Stop();
            Console.WriteLine("resolution done in {0} miliseconds", resTime.ElapsedMilliseconds);

            dmTime.Start();
            W.BuildDesignMatrix();
            dmTime.Stop();
            Console.WriteLine("design matrix done in {0} miliseconds", dmTime.ElapsedMilliseconds);

            abdTime.Start();
            W.FitAbundances();
            abdTime.Stop();
            Console.WriteLine("abundances done in {0} miliseconds", abdTime.ElapsedMilliseconds);

            Console.WriteLine("total time {0} seconds", (baselineTime.ElapsedMilliseconds + mofTime.ElapsedMilliseconds + resTime.ElapsedMilliseconds + dmTime.ElapsedMilliseconds + abdTime.ElapsedMilliseconds) / 1000d);

            //using (FileStream f = File.Open("abd.txt", FileMode.Create))
            //{
            //    using (StreamWriter sw = new StreamWriter(f))
            //    {
            //        for (int i = 0; i < W.Cluster.Count; i++)
            //        {
            //            sw.WriteLine(W.Cluster[i].CentreOfMass.ToString(bodka) + " " + W.Abundances[i].ToString(bodka));
            //        }                    
            //    } 
            //}

            //Vector<double> calcSpectrum = Vector<double>.Build.Dense((int)W.EndIndex);
            //double[] calcSpectrum = new double[(int)W.EndIndex];
            //Vector<double> abdVec = Vector<double>.Build.Dense(W.Abundances);

            //now we have to remove the last vector from the design matrix

            //double[] values = new double[W.DesignMatrix.Storage.ColumnPointers[W.DesignMatrix.Storage.ColumnPointers.Length - 2]];
            //int[] rowIndices = new int[W.DesignMatrix.Storage.ColumnPointers[W.DesignMatrix.Storage.ColumnPointers.Length - 2]];
            //int[] colPointers = new int[W.DesignMatrix.Storage.ColumnPointers.Length - 1];

            //Array.Copy(W.DesignMatrix.Storage.Values, values, values.Length);
            //Array.Copy(W.DesignMatrix.Storage.RowIndices, rowIndices, rowIndices.Length);
            //Array.Copy(W.DesignMatrix.Storage.ColumnPointers, colPointers, colPointers.Length);

            //CSparse.Double.SparseMatrix qwerty = new CSparse.Double.SparseMatrix(W.DesignMatrix.Storage.RowCount, W.DesignMatrix.Storage.ColumnCount - 1)
            //{
            //    Values = values,
            //    RowIndices = rowIndices,
            //    ColumnPointers = colPointers
            //};

            //qwerty.Multiply(W.Abundances, calcSpectrum);            

            //using (FileStream f = File.Open("calcSpectrum.txt", FileMode.Create))
            //{
            //    using (StreamWriter sw = new StreamWriter(f))
            //    {
            //        for (int i = 0; i < W.SpectralData.RawLength; i++)
            //        {
            //            sw.WriteLine(W.SpectralData.MassAxis[i].ToString(bodka) + " " + calcSpectrum[i].ToString(bodka) + " " + W.SpectralData.SignalAxis[i].ToString(bodka));
            //        }
            //    }
            //}

            Console.ReadKey();
        }

    }
}
