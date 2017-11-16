using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IsotopeFit.Numerics
{
    /// <summary>
    /// Base class that handles data interpolations and related data.
    /// </summary>
    public abstract class Interpolation
    {
        protected double[] xValues;
        protected double[] yValues;
        
        internal abstract double Evaluate(double x);
        
        internal abstract double[] Evaluate(double[] x);
        
        [Serializable]
        public class InterpolationException : Exception
        {
            public InterpolationException() { }
            public InterpolationException(string message) : base(message) { }
            public InterpolationException(string message, Exception inner) : base(message, inner) { }
            protected InterpolationException(
              System.Runtime.Serialization.SerializationInfo info,
              System.Runtime.Serialization.StreamingContext context) : base(info, context) { }
        }
    }
}
