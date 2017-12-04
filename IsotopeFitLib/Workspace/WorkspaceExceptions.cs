using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IsotopeFit
{

    [Serializable]
    public class WorkspaceException : Exception
    {
        public WorkspaceException() { }
        public WorkspaceException(string message) : base(message) { }
        public WorkspaceException(string message, Exception inner) : base(message, inner) { }
        protected WorkspaceException(
          System.Runtime.Serialization.SerializationInfo info,
          System.Runtime.Serialization.StreamingContext context) : base(info, context) { }
    }
}
