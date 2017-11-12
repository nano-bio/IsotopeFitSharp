using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IsotopeFit
{

    [Serializable]
    public class WorkspaceNotDefinedException : Exception
    {
        public WorkspaceNotDefinedException() { }
        public WorkspaceNotDefinedException(string message) : base(message) { }
        public WorkspaceNotDefinedException(string message, Exception inner) : base(message, inner) { }
        protected WorkspaceNotDefinedException(
          System.Runtime.Serialization.SerializationInfo info,
          System.Runtime.Serialization.StreamingContext context) : base(info, context) { }
    }
}
