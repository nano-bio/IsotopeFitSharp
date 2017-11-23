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

            Workspace W = new Workspace("testfiles\\testfile.ifd");
            time.Start();

            W.CorrectBaseline();

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

            Vector<double> calcSpectrum = Vector<double>.Build.Dense((int)W.EndIndex);
            Vector<double> abdVec = Vector<double>.Build.Dense(W.Abundances);

            //now we have to remove the last vector from the design matrix
            Matrix<double> mehehe = W.designMatrix.Storage.SubMatrix(0, W.designMatrix.Storage.RowCount, 0, W.designMatrix.Storage.ColumnCount - 1);

            calcSpectrum = mehehe.Multiply(abdVec);

            using (FileStream f = File.OpenWrite("calcSpectrum.txt"))
            {
                using (StreamWriter sw = new StreamWriter(f))
                {
                    for (int i = 0; i < W.SpectralData.RawLength; i++)
                    {
                        sw.WriteLine(W.SpectralData.MassOffsetCorrAxis[i].ToString() + " " + calcSpectrum[i].ToString());
                    }
                }
            }

            Console.WriteLine("done " + time.Elapsed);
            //Console.ReadKey();
        }

    }
}
