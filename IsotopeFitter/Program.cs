using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using IsotopeFit;

namespace IsotopeFitter
{
    class Program
    {
        //TODO: vymysliet, napisat, otestovat sposob rucneho zapisovania do datovych struktur

        static void Main(string[] args)
        {
            Workspace W = new Workspace();  // "testfile.ifd"

            W.Dummy();

            Console.ReadKey();
        }
    }
}
