﻿////////////////////////////////////////////////////////////////////////////////////////////////////
// file:	MathExtensions.cs
//
// summary:	Implements the mathematics extensions class
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace ActivationFnBenchmarks
{
    /// <summary>   The mathematics extensions. </summary>
    internal static class MathExtensions
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   A double extension method that converts a value to a float. </summary>
        ///
        /// <param name="value">    The value to act on. </param>
        ///
        /// <returns>   Value as a float. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        public static float ToFloat(this double value)
        {
            return (float)value;
        }
    }
}
