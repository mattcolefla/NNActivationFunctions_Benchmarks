////////////////////////////////////////////////////////////////////////////////////////////////////
// file:	ActivationFunctions.cs
//
// summary:	Implements the activation functions class
////////////////////////////////////////////////////////////////////////////////////////////////////

using System;
using System.Runtime.CompilerServices;
using BenchmarkDotNet.Attributes;
using Redzen.Numerics;


namespace ActivationFnBenchmarks
{
    /// <summary>   Activation Function Benchmarks. </summary>
    public sealed class ActivationFunctions
    {
        /// <summary>   The loops. </summary>
        const int Loops = 1000000;
        /// <summary>   The x coordinate. </summary>
        readonly double[] _x = new double[1000];
        /// <summary>   The f. </summary>
        readonly float[] _f = new float[1000];
        /// <summary>   The fudge d. </summary>
        const double FudgeD = 0.00001D;
        /// <summary>   The fudge f. </summary>
        const float FudgeF = 0.00001F;
        /// <summary>   threshold (left). </summary>
        const float tl = 0.001f;
        /// <summary>   threshold (right). </summary>
        const float tr = 0.999f;




        /// <summary>   Default constructor. </summary>
        public ActivationFunctions()
        {
            // Create some random Gaussian values as the inputs to the activation functions.
            var gaussian = new ZigguratGaussianSampler(0);
            for (var i = 0; i < _x.Length; i++)
            {
                _x[i] = gaussian.NextDouble(0, 2.0);
                _f[i] = (float)gaussian.NextDouble(0, 2.0);
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Logistic function steep double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double LogisticFunctionSteepDouble()
        {
            var a = 0.0;
            for (var i = 0; i < Loops; i++)
            {
                a = 1.0 / (1.0 + Math.Exp(-4.9 * _x[i % _x.Length]));
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Logistic approximant steep double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double LogisticApproximantSteepDouble()
        {
            var a = 0.0;
            for (var i = 0; i < Loops; i++)
            {
                a = 1.0 / (1.0 + Exp(-4.9 * _x[i % _x.Length]));
            }
            return a;
        }

        /// <summary>   Logistic sigmoid double. </summary>
        [Benchmark]
        public void LogisticSigmoidDouble()
        {
            var a = 0.0D;
            for (var i = 0; i < Loops; i++)
            {
                if (_x[i % _x.Length] < -40.0)
                    a = 0.0;
                if (_x[i % _x.Length] > 40.0)
                    a = 1.0;
                a = 1.0 / (1.0 + Math.Exp(-a));
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Soft sign double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double SoftSignDouble()
        {
            var a = 0.0;
            for (var i = 0; i < Loops; i++)
            {
                a = 0.5 + (_x[i % _x.Length] / (2.0 * (0.2 + Math.Abs(_x[i % _x.Length]))));
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Polynomial approximant double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double PolynomialApproximantDouble()
        {
            var a = 0.0;
            var x = 0.0;
            for (var i = 0; i < Loops; i++)
            {
                x = _x[i % _x.Length];
                x = x * 4.9;
                var x2 = x * x;
                var e = 1.0 + Math.Abs(x) + x2 * 0.555 + x2 * x2 * 0.143;
                var f = (x > 0) ? (1.0 / e) : e;
                a = 1.0 / (1.0 + f);
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Quadratic sigmoid double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double QuadraticSigmoidDouble()
        {
            var a = 0.0;
            const double t = 0.999;
            double y;

            for (var i = 0; i < Loops; i++)
            {
                double sign = Math.Sign(_x[i % _x.Length]);
                _x[i % _x.Length] = Math.Abs(_x[i % _x.Length]);

                if (_x[i % _x.Length] >= 0 && _x[i % _x.Length] < t)
                {
                    y = t - ((_x[i % _x.Length] - t) * (_x[i % _x.Length] - t));
                }
                else 
                {
                    y = t + (_x[i % _x.Length] - t) * FudgeD;
                }

                a = (y * sign * 0.5) + 0.5;
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Re lu double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double ReLuDouble()
        {
            var a = 0.0;
            double y;

            for (var i = 0; i < Loops; i++)
            {
                a = y = _x[i % _x.Length] > 0.0 ? _x[i % _x.Length] : 0.0;
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Leaky re lu double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double LeakyReLuDouble()
        {
            var a = 0.0;
            double y;

            for (var i = 0; i < Loops; i++)
            {
                if (_x[i % _x.Length] > 0.0)
                {
                    y = _x[i % _x.Length];
                }
                else
                {
                    y = _x[i % _x.Length] * FudgeD;
                }
                a = y;
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Leaky re lu shifted double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double LeakyReLuShiftedDouble()
        {
            var a = 0.0;
            for (var i = 0; i < Loops; i++)
            {
                const double offset = 0.5;

                double y;
                if (_x[i % _x.Length] + offset > 0.0)
                {
                    y = _x[i % _x.Length] + offset;
                }
                else
                {
                    y = (_x[i % _x.Length] + offset) * FudgeD;
                }
                a = y;
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Re lu double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double SreLuDouble()
        {
            var a = 0.0;
            const double tl = 0.001; // threshold (left).
            const double tr = 0.999; // threshold (right).

            for (var i = 0; i < Loops; i++)
            {
                double y;
                if (_x[i % _x.Length] > tl && _x[i % _x.Length] < tr)
                {
                    y = _x[i % _x.Length];
                }
                else if (_x[i % _x.Length] <= tl)
                {
                    y = tl + (_x[i % _x.Length] - tl) * a;
                }
                else
                {
                    y = tr + (_x[i % _x.Length] - tr) * a;
                }

                a = y;
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Re lu shifted double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double SreLuShiftedDouble()
        {
            var a = 0.0;
            const double tl = 0.001; // threshold (left).
            const double tr = 0.999; // threshold (right).
            const double offset = 0.5;

            for (var i = 0; i < Loops; i++)
            {
                double y;
                if (_x[i % _x.Length] + offset > tl && _x[i % _x.Length] + offset < tr)
                {
                    y = _x[i % _x.Length] + offset;
                }
                else if (_x[i % _x.Length] + offset <= tl)
                {
                    y = tl + ((_x[i % _x.Length] + offset) - tl) * a;
                }
                else
                {
                    y = tr + ((_x[i % _x.Length] + offset) - tr) * a;
                }

                a = y;
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Arc tangent double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double ArcTanDouble()
        {
            var a = 0.0;

            for (var i = 0; i < Loops; i++)
            {
                a = (Math.Atan(_x[i % _x.Length]) + Math.PI / 2.0) * 1.0 / Math.PI;
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Tangent h double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double TanHDouble()
        {
            var a = 0.0;
            for (var i = 0; i < Loops; i++)
            {
                a = (Math.Tanh(_x[i % _x.Length]) + 1.0) * 0.5;
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>
        /// Arc sine h double. Scaling factor from:
        /// https://www.reddit.com/r/MachineLearning/comments/6g5tg1/r_selfnormalizing_neural_networks_improved_elu/diwq7rb/.
        /// </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double ArcSinHDouble()
        {
            var a = 0.0;
            for (var i = 0; i < Loops; i++)
            {
                a = 1.2567348023993685 * ((Asinh(_x[i % _x.Length]) + 1.0) * 0.5);
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Scaled elu double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double ScaledEluDouble()
        {
            var a = 0.0;
            var alpha = 1.6732632423543772848170429916717;
            var scale = 1.0507009873554804934193349852946;

            for (var i = 0; i < Loops; i++)
            {
                double y;
                if (_x[i % _x.Length] >= 0)
                {
                    y = scale * _x[i % _x.Length];
                }
                else
                {
                    y = scale * ((alpha * Math.Exp(_x[i % _x.Length])) - alpha);
                }

                a = y;
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Maximum minus one double. </summary>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [Benchmark]
        public double MaxMinusOneDouble()
        {
            var a = 0.0;
            double y;

            for (var i = 0; i < Loops; i++)
            {
                if (_x[i % _x.Length] > -1)
                {
                    y = _x[i % _x.Length];
                }
                else
                {
                    y = -1;
                }
                a = y;
            }
            return a;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Exponent f. </summary>
        ///
        /// <param name="fVal"> The value. </param>
        ///
        /// <returns>   A float. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        private float ExpF(float fVal)
        {
            var tmp = (long) (1512775 * (double)fVal + (1072693248 - 60801));
            return (float) BitConverter.Int64BitsToDouble(tmp << 32);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Exps. </summary>
        ///
        /// <param name="val">  The value. </param>
        ///
        /// <returns>   A double. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Exp(double val)
        {
            long tmp = (long)(1512775 * val + (1072693248 - 60801));
            return BitConverter.Int64BitsToDouble(tmp << 32);
        }


        /// <summary>   Logistic function steep float. </summary>
        [Benchmark]
        public void LogisticFunctionSteepFloat()
        {
            var f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                f = 1.0f / (1.0f + (float)ExpF(-4.9f * _f[i % _f.Length]));
            }
        }

        /// <summary>   Logistic approximant steep float. </summary>
        [Benchmark]
        public void LogisticApproximantSteepFloat()
        {
            var f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                f = 1f / (1f + ExpF(-4.9f * _f[i % _f.Length]));
            }
        }

        /// <summary>   Logistic sigmoid float. </summary>
        [Benchmark]
        public void LogisticSigmoidFloat()
        {
            var f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                if (_f[i % _f.Length] < -40.0F)
                    f = 0.0F;
                if (_f[i % _f.Length] > 40.0F)
                    f = 1.0F;
                f = 1.0F / (1.0F + (float)ExpF(-_f[i % _f.Length]));
            }
        }

        /// <summary>   Soft sign float. </summary>
        [Benchmark]
        public void SoftSignFloat()
        {
            var f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                f = 0.5f + (_f[i % _f.Length] / (2.0f * (0.2f + Math.Abs(_f[i % _f.Length]))));
            }
        }

        /// <summary>   Polynomial approximant float. </summary>
        [Benchmark]
        public void PolynomialApproximantFloat()
        {
            var x2 = 0.0F;
            var e = 0.0F;
            var f = 0.0F;
            var x = 0.0F;

            for (var i = 0; i < Loops; i++)
            {
                x = _f[i % _f.Length];

                // Very close approximation to LogisticFunctionSteep that avoids exp.
                x = x * 4.9f;
                x2 = x * x;
                e = 1.0f + Math.Abs(x) + x2 * 0.555f + x2 * x2 * 0.143f;
                f = (x > 0f) ? (1.0f / e) : e;
                f =1.0f / (1.0f + f);
            }
        }

        /// <summary>   Quadratic sigmoid float. </summary>
        [Benchmark]
        public void QuadraticSigmoidFloat()
        {
            const float t = 0.999F;
            var f = 0.0F;

            for (var i = 0; i < Loops; i++)
            {
                f = _f[i % _f.Length];

                float sign = Math.Sign(_f[i % _f.Length]);
                _f[i % _f.Length] = Math.Abs(_f[i % _f.Length]);

                float y = 0;
                if (_f[i % _f.Length] >= 0 && _f[i % _f.Length] < t)
                {
                    y = t - ((_f[i % _f.Length] - t) * (_f[i % _f.Length] - t));
                }
                else //if (x >= t) 
                {
                    y = t + (_f[i % _f.Length] - t) * FudgeF;
                }

                f = (y * sign * 0.5f) + 0.5f;
            }
        }

        /// <summary>   Re lu float. </summary>
        [Benchmark]
        public void SreLuFloat()
        {
            var f = 0.0F;

            for (var i = 0; i < Loops; i++)
            {
                float y;
                if (_f[i % _f.Length] > tl && _f[i % _f.Length] < tr)
                {
                    y = _f[i % _f.Length];
                }
                else if (_f[i % _f.Length] <= tl)
                {
                    y = tl + (_f[i % _f.Length] - tl) * FudgeF;
                }
                else
                {
                    y = tr + (_f[i % _f.Length] - tr) * FudgeF;
                }

                f = y;
            }
        }

        /// <summary>   Re lu shifted float. </summary>
        [Benchmark]
        public void SreLuShiftedFloat()
        {
            float f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                if (i == 0)
                    f = _f[i % _f.Length] += 0.5f;

                float y;
                if (_f[i % _f.Length] > tl && _f[i % _f.Length] < tr)
                {
                    y = _f[i % _f.Length];
                }
                else if (_f[i % _f.Length] <= tl)
                {
                    y = tl + (_f[i % _f.Length] - tl) * FudgeF;
                }
                else
                {
                    y = tr + (_f[i % _f.Length] - tr) * FudgeF;
                }

                f = y;
            }
        }

        /// <summary>   Re lu float. </summary>
        [Benchmark]
        public void ReLuFloat()
        {
            float y = 0.0F;

            for (var i = 0; i < Loops; i++)
            {
                y = _f[i % _f.Length] > 0.0F ? _f[i % _f.Length] : 0.0F;
            }
        }

        /// <summary>   Leaky re lu float. </summary>
        [Benchmark]
        public void LeakyReLuFloat()
        {
            for (var i = 0; i < Loops; i++)
            {
                float y;
                if (_f[i % _f.Length] > 0.0F)
                {
                    y = _f[i % _f.Length];
                }
                else
                {
                    y = _f[i % _f.Length] * FudgeF;
                }
            }
        }

        /// <summary>   Leaky re lu shifted float. </summary>
        [Benchmark]
        public void LeakyReLuShiftedFloat()
        {
            const float offset = 0.5F;

            for (var i = 0; i < Loops; i++)
            {
                float y = (_f[i % _f.Length] + offset > 0.0F) ? _f[i % _f.Length] + offset : (_f[i % _f.Length] + offset) * FudgeF;
            }
        }


        /// <summary>   Arc tangent float. </summary>
        [Benchmark]
        public void ArcTanFloat()
        {
            float f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                f =(Math.Atan(_f[i % _f.Length]).ToFloat() + (float)Math.PI / 2.0F) * 1.0F / (float)Math.PI;
            }
        }

        /// <summary>   Tangent h float. </summary>
        [Benchmark]
        public void TanHFloat()
        {
            float f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                // Scaling factor from:
                // https://www.reddit.com/r/MachineLearning/comments/6g5tg1/r_selfnormalizing_neural_networks_improved_elu/diwq7rb/

                f = (Math.Tanh(_f[i % _f.Length]).ToFloat() + 1.0F) * 0.5F;
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Asinh f. </summary>
        ///
        /// <param name="value">    The value. </param>
        ///
        /// <returns>   A float. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float AsinhF(float value)
        {
            FastLog.Log2(1); //create table for the fist time
            float mathEf = 2.7182818284590451F;
            return FastLog.Log(value + Math.Sqrt((value * value) + 1F).ToFloat(), mathEf);
        }

        /// <summary>   Arc sine h float. </summary>
        [Benchmark]
        public void ArcSinHFloat()
        {
            var f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                // Scaling factor from:
                // https://www.reddit.com/r/MachineLearning/comments/6g5tg1/r_selfnormalizing_neural_networks_improved_elu/diwq7rb/
                f = 1.2567348023993685F * ((AsinhF(_f[i % _f.Length]) + 1.0F) * 0.5F);
            }
        }

        /// <summary>   Scaled elu float. </summary>
        [Benchmark]
        public void ScaledEluFloat()
        {
            float alpha = 1.6732632423543772848170429916717F;
            float scale = 1.0507009873554804934193349852946F;
            float f = 0.0F;

            for (var i = 0; i < Loops; i++)
            {
                float y;
                if (_f[i % _f.Length] >= 0F)
                {
                    y = scale * _f[i % _f.Length];
                }
                else
                {
                    y = scale * ((alpha * ExpF(_f[i % _f.Length])) - alpha);
                }
            }
        }

        /// <summary>   Maximum minus one float. </summary>
        [Benchmark]
        public void MaxMinusOneFloat()
        {
            for (var i = 0; i < Loops; i++)
            {
                float y;
                if (_f[i % _f.Length] > -1F)
                {
                    y = _f[i % _f.Length];
                }
                else
                {
                    y = -1F;
                }
            }
        }

        /// <summary>   new. </summary>
        [Benchmark]
        public void BinaryStepFloat()
        {
            float f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                f = _f[i % _f.Length] < 0F ? 0F : 1F;
            }
        }

        /// <summary>   Binary step double. </summary>
        [Benchmark]
        public void BinaryStepDouble()
        {
            float f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                f = _f[i % _f.Length] < 0 ? 0 : 1;
            }
        }


        /// <summary>   Parameteric re lu float. </summary>
        [Benchmark]
        public void ParametericReLuFloat()
        {
            float f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                f = _f[i % _f.Length] < 0F ? FudgeF * _f[i % _f.Length] : _f[i % _f.Length];
            }
        }

        /// <summary>   Parameteric re lu double. </summary>
        [Benchmark]
        public void ParametericReLuDouble()
        {
            double f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                f = _x[i % _x.Length] < 0 ? 5 * _x[i % _x.Length] : _x[i % _x.Length];
            }
        }


        /// <summary>   Bent identity float. </summary>
        [Benchmark]
        public void BentIdentityFloat()
        {
            float f = 0.0F;

            for (var i = 0; i < Loops; i++)
            {
                f = (float) (((Math.Sqrt(Math.Pow(_f[i % _f.Length], 2F) + 1F) - 1F) / 2F) + _f[i % _f.Length]);
            }
        }

        /// <summary>   Bent identity double. </summary>
        [Benchmark]
        public void BentIdentityDouble()
        {
            double f = 0.0F;
            for (var i = 0; i < Loops; i++)
            {
                f = (((Math.Sqrt(Math.Pow(_x[i % _x.Length], 2) + 1)) - 1) / 2) + _x[i % _x.Length];
            }
        }



        #region Private Static Methods

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /// <summary>   Hyperbolic Area Sine. </summary>
        ///
        /// <param name="value">    The real value. </param>
        ///
        /// <returns>   The hyperbolic angle, i.e. the area of its hyperbolic sector. </returns>
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Asinh(double value)
        {
            return Math.Log(value + Math.Sqrt((value * value) + 1), Math.E);
        }

        #endregion
    }
}
