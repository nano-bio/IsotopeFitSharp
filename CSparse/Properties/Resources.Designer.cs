﻿//------------------------------------------------------------------------------
// <auto-generated>
//     This code was generated by a tool.
//     Runtime Version:4.0.30319.34209
//
//     Changes to this file may cause incorrect behavior and will be lost if
//     the code is regenerated.
// </auto-generated>
//------------------------------------------------------------------------------

namespace CSparse.Properties {
    using System;
    
    
    /// <summary>
    ///   A strongly-typed resource class, for looking up localized strings, etc.
    /// </summary>
    // This class was auto-generated by the StronglyTypedResourceBuilder
    // class via a tool like ResGen or Visual Studio.
    // To add or remove a member, edit your .ResX file then rerun ResGen
    // with the /str option, or rebuild your VS project.
    [global::System.CodeDom.Compiler.GeneratedCodeAttribute("System.Resources.Tools.StronglyTypedResourceBuilder", "4.0.0.0")]
    [global::System.Diagnostics.DebuggerNonUserCodeAttribute()]
    [global::System.Runtime.CompilerServices.CompilerGeneratedAttribute()]
    internal class Resources {
        
        private static global::System.Resources.ResourceManager resourceMan;
        
        private static global::System.Globalization.CultureInfo resourceCulture;
        
        [global::System.Diagnostics.CodeAnalysis.SuppressMessageAttribute("Microsoft.Performance", "CA1811:AvoidUncalledPrivateCode")]
        internal Resources() {
        }
        
        /// <summary>
        ///   Returns the cached ResourceManager instance used by this class.
        /// </summary>
        [global::System.ComponentModel.EditorBrowsableAttribute(global::System.ComponentModel.EditorBrowsableState.Advanced)]
        internal static global::System.Resources.ResourceManager ResourceManager {
            get {
                if (object.ReferenceEquals(resourceMan, null)) {
                    global::System.Resources.ResourceManager temp = new global::System.Resources.ResourceManager("CSparse.Properties.Resources", typeof(Resources).Assembly);
                    resourceMan = temp;
                }
                return resourceMan;
            }
        }
        
        /// <summary>
        ///   Overrides the current thread's CurrentUICulture property for all
        ///   resource lookups using this strongly typed resource class.
        /// </summary>
        [global::System.ComponentModel.EditorBrowsableAttribute(global::System.ComponentModel.EditorBrowsableState.Advanced)]
        internal static global::System.Globalization.CultureInfo Culture {
            get {
                return resourceCulture;
            }
            set {
                resourceCulture = value;
            }
        }
        
        /// <summary>
        ///   Looks up a localized string similar to Invalid column ordering..
        /// </summary>
        internal static string InvalidColumnOrdering {
            get {
                return ResourceManager.GetString("InvalidColumnOrdering", resourceCulture);
            }
        }
        
        /// <summary>
        ///   Looks up a localized string similar to Invalid matrix dimensions..
        /// </summary>
        internal static string InvalidDimensions {
            get {
                return ResourceManager.GetString("InvalidDimensions", resourceCulture);
            }
        }
        
        /// <summary>
        ///   Looks up a localized string similar to Invalid permutation vector..
        /// </summary>
        internal static string InvalidPermutation {
            get {
                return ResourceManager.GetString("InvalidPermutation", resourceCulture);
            }
        }
        
        /// <summary>
        ///   Looks up a localized string similar to Matrix dimension must not be a negative number..
        /// </summary>
        internal static string MatrixDimensionNonNegative {
            get {
                return ResourceManager.GetString("MatrixDimensionNonNegative", resourceCulture);
            }
        }
        
        /// <summary>
        ///   Looks up a localized string similar to Matrix dimension must be a positive number..
        /// </summary>
        internal static string MatrixDimensionPositive {
            get {
                return ResourceManager.GetString("MatrixDimensionPositive", resourceCulture);
            }
        }
        
        /// <summary>
        ///   Looks up a localized string similar to Matrix dimensions don&apos;t match..
        /// </summary>
        internal static string MatrixDimensions {
            get {
                return ResourceManager.GetString("MatrixDimensions", resourceCulture);
            }
        }
        
        /// <summary>
        ///   Looks up a localized string similar to Matrix must be square..
        /// </summary>
        internal static string MatrixSquare {
            get {
                return ResourceManager.GetString("MatrixSquare", resourceCulture);
            }
        }
        
        /// <summary>
        ///   Looks up a localized string similar to Matrix must be symmetric positive definite..
        /// </summary>
        internal static string MatrixSymmetricPositiveDefinite {
            get {
                return ResourceManager.GetString("MatrixSymmetricPositiveDefinite", resourceCulture);
            }
        }
        
        /// <summary>
        ///   Looks up a localized string similar to Value must not be NaN..
        /// </summary>
        internal static string ValueNotNaN {
            get {
                return ResourceManager.GetString("ValueNotNaN", resourceCulture);
            }
        }
        
        /// <summary>
        ///   Looks up a localized string similar to Vectors must have the same dimension..
        /// </summary>
        internal static string VectorsSameLength {
            get {
                return ResourceManager.GetString("VectorsSameLength", resourceCulture);
            }
        }
    }
}
