using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.Diagnostics;

using IsotopeFit;

namespace IsotopeFitter
{
    class Program
    {
        //TODO: vymysliet, napisat, otestovat sposob rucneho zapisovania do datovych struktur

        static void Main(string[] args)
        {
            Workspace W = new Workspace("testfile.ifd");  // 

            int idx = 0;

            List<double> axisList = W.RawData.MassAxis.ToList();
            axisList.Sort();

            double[] axis = axisList.ToArray();

            Stopwatch sw = new Stopwatch();

            sw.Start();
            idx = FindLowerLimitIndex(axis, 3000);
            sw.Stop();

            double elapsed = sw.ElapsedTicks / Stopwatch.Frequency;

            Console.WriteLine(elapsed);
            Console.WriteLine(idx + " " + axis[idx]);

            Console.ReadKey();
        }

        static int FindLowerLimitIndex(double[] array, double threshold)
        {
            int index = Array.BinarySearch(array, threshold);

            if (index >= 0) return index;
            else return ~index;
        }
    }
}
