using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Reflection;

using IsotopeFit;

namespace IsotopeFitLib.Tests
{
    [TestFixture]
    public partial class Tests
    {
        [Test]
        public void DesignMatrixBuildTest()
        {
            Workspace wrk = new Workspace(Path.GetDirectoryName(Assembly.GetAssembly(typeof(Tests)).Location) + "\\TestData\\testfile.ifd");

            wrk.BuildDesignMatrix();

            Assert.Warn("Design matrix build test passed, but this does not guarantee that the design matrix is correct!");
            //Assert.Pass("Design matrix build test passed.");
        }
    }
}
