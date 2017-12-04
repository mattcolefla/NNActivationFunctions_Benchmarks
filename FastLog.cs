////////////////////////////////////////////////////////////////////////////////////////////////////
// file:	FastLog.cs
//
// summary:	Implements the fast log class
////////////////////////////////////////////////////////////////////////////////////////////////////

using System;
using System.Runtime.InteropServices;

namespace ActivationFnBenchmarks
{
    /// <summary>   A fast log. </summary>
    internal static class FastLog
    {
        /// <summary>   An ieee 754. </summary>
        [StructLayout(LayoutKind.Explicit)]
        private struct Ieee754
        {
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            /// <summary>   Gets the single. </summary>
            ///
            /// <value> The single. </value>
            ////////////////////////////////////////////////////////////////////////////////////////////////////

            [FieldOffset(0)] public float Single;

            ////////////////////////////////////////////////////////////////////////////////////////////////////
            /// <summary>   Gets the unsigned bits. </summary>
            ///
            /// <value> The unsigned bits. </value>
            ////////////////////////////////////////////////////////////////////////////////////////////////////

            [FieldOffset(0)] public uint UnsignedBits;

            ////////////////////////////////////////////////////////////////////////////////////////////////////
            /// <summary>   Gets the signed bits. </summary>
            ///
            /// <value> The signed bits. </value>
            ////////////////////////////////////////////////////////////////////////////////////////////////////

            [FieldOffset(0)] public int SignedBits;

            ////////////////////////////////////////////////////////////////////////////////////////////////////
            /// <summary>   Gets the sign. </summary>
            ///
            /// <value> The sign. </value>
            ////////////////////////////////////////////////////////////////////////////////////////////////////

            public uint Sign => UnsignedBits >> 31;

            ////////////////////////////////////////////////////////////////////////////////////////////////////
            /// <summary>   Gets the exponent. </summary>
            ///
            /// <value> The exponent. </value>
            ////////////////////////////////////////////////////////////////////////////////////////////////////

            public int Exponent => (SignedBits >> 23) & 0xFF;

            ////////////////////////////////////////////////////////////////////////////////////////////////////
            /// <summary>   Gets the mantissa. </summary>
            ///
            /// <value> The mantissa. </value>
            ////////////////////////////////////////////////////////////////////////////////////////////////////

            public uint Mantissa => UnsignedBits & 0x007FFFFF;
        }

        /// <summary>   The mantissa logs. </summary>
        private static readonly float[] MantissaLogs = new float[(int)Math.Pow(2, 23)];
        /// <summary>   The base 10. </summary>
        private const float Base10 = 3.321928F;
        /// <summary>   The base. </summary>
        private const float BaseE = 1.442695F;

        /// <summary>   Static constructor. </summary>
        static FastLog()
        {
            //creating lookup table
            for (uint i = 0; i < MantissaLogs.Length; i++)
            {
                var n = new Ieee754 { UnsignedBits = i | 0x3F800000 }; //added the implicit 1 leading bit
                MantissaLogs[i] = (float)Math.Log(n.Single, 2);
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Logs a 2. </summary>
        ///
        /// <param name="value">    The value. </param>
        ///
        /// <returns>   A float. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        public static float Log2(float value)
        {
            if (Math.Abs(value) < 0F)
            {
                return float.NegativeInfinity;
            }

            var number = new Ieee754 { Single = value };

            if (number.UnsignedBits >> 31 == 1) //NOTE: didn't call Sign property for higher performance
                return float.NaN;

            return (((number.SignedBits >> 23) & 0xFF) - 127) + MantissaLogs[number.UnsignedBits & 0x007FFFFF];
            //NOTE: didn't call Exponent and Mantissa properties for higher performance
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Logs a 10. </summary>
        ///
        /// <param name="value">    The value. </param>
        ///
        /// <returns>   A float. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        public static float Log10(float value)
        {
            return Log2(value) / Base10;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Lns. </summary>
        ///
        /// <param name="value">    The value. </param>
        ///
        /// <returns>   A float. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        public static float Ln(float value)
        {
            return Log2(value) / BaseE;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Logs. </summary>
        ///
        /// <param name="value">        The value. </param>
        /// <param name="valueBase">    The value base. </param>
        ///
        /// <returns>   A float. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        public static float Log(float value, float valueBase)
        {
            return Log2(value) / Log2(valueBase);
        }
    }
}
