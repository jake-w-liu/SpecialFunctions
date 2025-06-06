using System;
using System.Collections.Generic;
using System.Globalization;
using System.Runtime.Serialization;
using System.Linq;
using Complex = System.Numerics.Complex;
using System.Text;
#if NET5_0_OR_GREATER
using System.Text.Json.Serialization;
#endif

using System.Runtime;

namespace SpecialFunctions
{

    [Serializable]
    [DataContract(Namespace = "urn:MathNet/Numerics")]
    public class Polynomial : IFormattable, IEquatable<Polynomial>, ICloneable
    {

        [DataMember(Order = 1)]
#if NET5_0_OR_GREATER
        [JsonInclude]
#endif
        public double[] Coefficients { get; private set; }


        [DataMember(Order = 2)]
        public string VariableName = "x";

        public int Degree => EvaluateDegree(Coefficients);

        public Polynomial(int n)
        {
            if (n < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(n), "n must be non-negative");
            }

            Coefficients = new double[n];
        }

        /// <summary>
        /// Create a zero-polynomial
        /// </summary>
        public Polynomial()
        {
            Coefficients = Array.Empty<double>();
        }

        public Polynomial(double coefficient)
        {
            if (coefficient == 0.0)
            {
                Coefficients = Array.Empty<double>();
            }
            else
            {
                Coefficients = new[] { coefficient };
            }
        }

        public Polynomial(params double[] coefficients)
        {
            Coefficients = coefficients;
        }

        public Polynomial(IEnumerable<double> coefficients) : this(coefficients.ToArray())
        {
        }

        public static Polynomial Zero => new Polynomial();


        static int EvaluateDegree(double[] coefficients)
        {
            for (int i = coefficients.Length - 1; i >= 0; i--)
            {
                if (coefficients[i] != 0.0)
                {
                    return i;
                }
            }

            return -1;
        }

        #region Evaluation

        public static double Evaluate(double z, params double[] coefficients)
        {

            // 2020-10-07 jbialogrodzki #730 Since this is public API we should probably
            // handle null arguments? It doesn't seem to have been done consistently in this class though.
            if (coefficients == null)
            {
                throw new ArgumentNullException(nameof(coefficients));
            }

            // 2020-10-07 jbialogrodzki #730 Zero polynomials need explicit handling.
            // Without this check, we attempted to peek coefficients at negative indices!
            int n = coefficients.Length;
            if (n == 0)
            {
                return 0;
            }

            double sum = coefficients[n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                sum *= z;
                sum += coefficients[i];
            }

            return sum;

        }

     
        public static Complex Evaluate(Complex z, params double[] coefficients)
        {

            // 2020-10-07 jbialogrodzki #730 Since this is a public API we should probably
            // handle null arguments? It doesn't seem to have been done consistently in this class though.
            if (coefficients == null)
            {
                throw new ArgumentNullException(nameof(coefficients));
            }

            // 2020-10-07 jbialogrodzki #730 Zero polynomials need explicit handling.
            // Without this check, we attempted to peek coefficients at negative indices!
            int n = coefficients.Length;
            if (n == 0)
            {
                return 0;
            }

            Complex sum = coefficients[n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                sum *= z;
                sum += coefficients[i];
            }

            return sum;

        }

  
        public static Complex Evaluate(Complex z, params Complex[] coefficients)
        {

            // 2020-10-07 jbialogrodzki #730 Since this is a public API we should probably
            // handle null arguments? It doesn't seem to have been done consistently in this class though.
            if (coefficients == null)
            {
                throw new ArgumentNullException(nameof(coefficients));
            }

            // 2020-10-07 jbialogrodzki #730 Zero polynomials need explicit handling.
            // Without this check, we attempted to peek coefficients at negative indices!
            int n = coefficients.Length;
            if (n == 0)
            {
                return 0;
            }

            Complex sum = coefficients[n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                sum *= z;
                sum += coefficients[i];
            }

            return sum;

        }

        public double Evaluate(double z)
        {
            return Evaluate(z, Coefficients);
        }

        public Complex Evaluate(Complex z)
        {
            return Evaluate(z, Coefficients);
        }

        /// <summary>
        /// Evaluate a polynomial at points z.
        /// </summary>
        /// <param name="z">The locations where to evaluate the polynomial at.</param>
        public IEnumerable<double> Evaluate(IEnumerable<double> z)
        {
            return z.Select(Evaluate);
        }

        public IEnumerable<Complex> Evaluate(IEnumerable<Complex> z)
        {
            return z.Select(Evaluate);
        }

        #endregion

        #region Calculus

        public Polynomial Differentiate()
        {
            int n = Degree;
            if (n < 0)
            {
                return this;
            }

            if (n == 0)
            {
                // Zero
                return Zero;
            }

            var c = new double[n];
            for (int i = 0; i < c.Length; i++)
            {
                c[i] = Coefficients[i + 1] * (i + 1);
            }

            return new Polynomial(c);
        }

        public Polynomial Integrate()
        {
            int n = Degree;
            if (n < 0)
            {
                return this;
            }

            var c = new double[n + 2];
            for (int i = 1; i < c.Length; i++)
            {
                c[i] = Coefficients[i - 1] / i;
            }

            return new Polynomial(c);
        }

        #endregion


        #region Arithmetic Operations

        public static Polynomial Add(Polynomial a, Polynomial b)
        {
            var ac = a.Coefficients;
            var bc = b.Coefficients;

            var degree = Math.Max(a.Degree, b.Degree);
            var result = new double[degree + 1];

            var commonLength = Math.Min(Math.Min(ac.Length, bc.Length), result.Length);
            for (int i = 0; i < commonLength; i++)
            {
                result[i] = ac[i] + bc[i];
            }

            int acLength = Math.Min(ac.Length, result.Length);
            for (int i = commonLength; i < acLength; i++)
            {
                // no need to add since only one of both applies
                result[i] = ac[i];
            }

            int bcLength = Math.Min(bc.Length, result.Length);
            for (int i = commonLength; i < bcLength; i++)
            {
                // no need to add since only one of both applies
                result[i] = bc[i];
            }

            return new Polynomial(result);
        }

        public static Polynomial Add(Polynomial a, double b)
        {
            var ac = a.Coefficients;

            var degree = Math.Max(a.Degree, 0);
            var result = new double[degree + 1];

            var commonLength = Math.Min(ac.Length, result.Length);
            for (int i = 0; i < commonLength; i++)
            {
                result[i] = ac[i];
            }

            result[0] += b;

            return new Polynomial(result);
        }

        public static Polynomial Subtract(Polynomial a, Polynomial b)
        {
            var ac = a.Coefficients;
            var bc = b.Coefficients;

            var degree = Math.Max(a.Degree, b.Degree);
            var result = new double[degree + 1];

            var commonLength = Math.Min(Math.Min(ac.Length, bc.Length), result.Length);
            for (int i = 0; i < commonLength; i++)
            {
                result[i] = ac[i] - bc[i];
            }

            int acLength = Math.Min(ac.Length, result.Length);
            for (int i = commonLength; i < acLength; i++)
            {
                // no need to add since only one of both applies
                result[i] = ac[i];
            }

            int bcLength = Math.Min(bc.Length, result.Length);
            for (int i = commonLength; i < bcLength; i++)
            {
                // no need to add since only one of both applies
                result[i] = -bc[i];
            }

            return new Polynomial(result);
        }
        public static Polynomial Subtract(Polynomial a, double b)
        {
            return Add(a, -b);
        }
        public static Polynomial Subtract(double b, Polynomial a)
        {
            var ac = a.Coefficients;

            var degree = Math.Max(a.Degree, 0);
            var result = new double[degree + 1];

            var commonLength = Math.Min(ac.Length, result.Length);
            for (int i = 0; i < commonLength; i++)
            {
                result[i] = -ac[i];
            }

            result[0] += b;

            return new Polynomial(result);
        }

  
        public static Polynomial Negate(Polynomial a)
        {
            var ac = a.Coefficients;

            var degree = a.Degree;
            var result = new double[degree + 1];

            for (int i = 0; i < result.Length; i++)
            {
                result[i] = -ac[i];
            }

            return new Polynomial(result);
        }

        public static Polynomial Multiply(Polynomial a, Polynomial b)
        {

            // 2020-10-07 jbialogrodzki #730 Since this is a public API we should probably
            // handle null arguments? It doesn't seem to have been done consistently in this class though.
            if (a == null)
            {
                throw new ArgumentNullException(nameof(a));
            }

            if (b == null)
            {
                throw new ArgumentNullException(nameof(b));
            }

            var ad = a.Degree;
            var bd = b.Degree;

            // 2020-10-07 jbialogrodzki #730 Zero polynomials need explicit handling.
            // Without this check, we attempted to create arrays of negative lengths!
            if (ad < 0 || bd < 0)
            {
                return Polynomial.Zero;
            }

            double[] ac = a.Coefficients;
            double[] bc = b.Coefficients;

            var degree = ad + bd;
            double[] result = new double[degree + 1];

            for (int i = 0; i <= ad; i++)
            {
                for (int j = 0; j <= bd; j++)
                {
                    result[i + j] += ac[i] * bc[j];
                }
            }

            return new Polynomial(result);

        }

        public static Polynomial Multiply(Polynomial a, double k)
        {
            var ac = a.Coefficients;

            var result = new double[a.Degree + 1];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = ac[i] * k;
            }

            return new Polynomial(result);
        }

      
        public static Polynomial Divide(Polynomial a, double k)
        {
            var ac = a.Coefficients;

            var result = new double[a.Degree + 1];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = ac[i] / k;
            }

            return new Polynomial(result);
        }

        public static (Polynomial, Polynomial) DivideRemainder(Polynomial a, Polynomial b)
        {
            var bDegree = b.Degree;
            if (bDegree < 0)
            {
                throw new DivideByZeroException("b polynomial ends with zero");
            }

            var aDegree = a.Degree;
            if (aDegree < 0)
            {
                // zero divided by non-zero is zero without remainder
                return (a, a);
            }

            if (bDegree == 0)
            {
                // division by scalar
                return (Divide(a, b.Coefficients[0]), Zero);
            }

            if (aDegree < bDegree)
            {
                // denominator degree higher than nominator degree
                // quotient always be 0 and return c1 as remainder
                return (Zero, a);
            }

            var c1 = a.Coefficients.ToArray();
            var c2 = b.Coefficients.ToArray();

            var scl = c2[bDegree];
            var c22 = new double[bDegree];
            for (int ii = 0; ii < c22.Length; ii++)
            {
                c22[ii] = c2[ii] / scl;
            }

            int i = aDegree - bDegree;
            int j = aDegree;
            while (i >= 0)
            {
                var v = c1[j];
                for (int k = i; k < j; k++)
                {
                    c1[k] -= c22[k - i] * v;
                }

                i--;
                j--;
            }

            var j1 = j + 1;
            var l1 = aDegree - j;

            var quo = new double[l1];
            for (int k = 0; k < l1; k++)
            {
                quo[k] = c1[k + j1] / scl;
            }

            var rem = new double[j1];
            for (int k = 0; k < j1; k++)
            {
                rem[k] = c1[k];
            }

            return (new Polynomial(quo), new Polynomial(rem));
        }

        #endregion

        #region Arithmetic Pointwise Operations

        public static Polynomial PointwiseDivide(Polynomial a, Polynomial b)
        {
            var ac = a.Coefficients;
            var bc = b.Coefficients;

            var degree = a.Degree;
            var result = new double[degree + 1];

            var commonLength = Math.Min(Math.Min(ac.Length, bc.Length), result.Length);
            for (int i = 0; i < commonLength; i++)
            {
                result[i] = ac[i] / bc[i];
            }

            for (int i = commonLength; i < result.Length; i++)
            {
                result[i] = ac[i] / 0.0;
            }

            return new Polynomial(result);
        }

        public static Polynomial PointwiseMultiply(Polynomial a, Polynomial b)
        {
            var ac = a.Coefficients;
            var bc = b.Coefficients;

            var degree = Math.Min(a.Degree, b.Degree);
            var result = new double[degree + 1];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = ac[i] * bc[i];
            }

            return new Polynomial(result);
        }

        #endregion

        #region Arithmetic Instance Methods (forwarders)

        public (Polynomial, Polynomial) DivideRemainder(Polynomial b)
        {
            return DivideRemainder(this, b);
        }

        #endregion

        #region Arithmetic Operator Overloads (forwarders)


        public static Polynomial operator +(Polynomial a, Polynomial b)
        {
            return Add(a, b);
        }


        public static Polynomial operator +(Polynomial a, double k)
        {
            return Add(a, k);
        }


        public static Polynomial operator +(double k, Polynomial a)
        {
            return Add(a, k);
        }


        public static Polynomial operator -(Polynomial a, Polynomial b)
        {
            return Subtract(a, b);
        }


        public static Polynomial operator -(Polynomial a, double k)
        {
            return Subtract(a, k);
        }

        public static Polynomial operator -(double k, Polynomial a)
        {
            return Subtract(k, a);
        }

        public static Polynomial operator -(Polynomial a)
        {
            return Negate(a);
        }

        public static Polynomial operator *(Polynomial a, Polynomial b)
        {
            return Multiply(a, b);
        }

        public static Polynomial operator *(Polynomial a, double k)
        {
            return Multiply(a, k);
        }

        public static Polynomial operator *(double k, Polynomial a)
        {
            return Multiply(a, k);
        }

        public static Polynomial operator /(Polynomial a, double k)
        {
            return Divide(a, k);
        }

        #endregion

        #region ToString

        public override string ToString()
        {
            return ToString("G", CultureInfo.CurrentCulture);
        }

        public string ToStringDescending()
        {
            return ToStringDescending("G", CultureInfo.CurrentCulture);
        }

        public string ToString(string format)
        {
            return ToString(format, CultureInfo.CurrentCulture);
        }
        public string ToStringDescending(string format)
        {
            return ToStringDescending(format, CultureInfo.CurrentCulture);
        }

        public string ToString(IFormatProvider formatProvider)
        {
            return ToString("G", formatProvider);
        }

        public string ToStringDescending(IFormatProvider formatProvider)
        {
            return ToStringDescending("G", formatProvider);
        }

#pragma warning disable CS8767 // Nullability of reference types in type of parameter doesn't match implicitly implemented member (possibly because of nullability attributes).
        public string ToString(string format, IFormatProvider formatProvider)
#pragma warning restore CS8767 // Nullability of reference types in type of parameter doesn't match implicitly implemented member (possibly because of nullability attributes).
        {
            if (Degree < 0)
            {
                return "0";
            }

            var sb = new StringBuilder();
            bool first = true;
            for (int i = 0; i < Coefficients.Length; i++)
            {
                double c = Coefficients[i];
                if (c == 0.0)
                {
                    continue;
                }

                if (first)
                {
                    sb.Append(c.ToString(format, formatProvider));
                    if (i > 0)
                    {
                        sb.Append(VariableName);
                    }

                    if (i > 1)
                    {
                        sb.Append("^");
                        sb.Append(i);
                    }

                    first = false;
                }
                else
                {
                    if (c < 0.0)
                    {
                        sb.Append(" - ");
                        sb.Append((-c).ToString(format, formatProvider));
                    }
                    else
                    {
                        sb.Append(" + ");
                        sb.Append(c.ToString(format, formatProvider));
                    }

                    if (i > 0)
                    {
                        sb.Append(VariableName);
                    }

                    if (i > 1)
                    {
                        sb.Append("^");
                        sb.Append(i);
                    }
                }
            }

            return sb.ToString();
        }

        /// <summary>
        /// Format the polynomial in descending order, e.g. "x^3 + 2.0x^2 - 4.3".
        /// </summary>
        public string ToStringDescending(string format, IFormatProvider formatProvider)
        {
            if (Degree < 0)
            {
                return "0";
            }

            var sb = new StringBuilder();
            bool first = true;
            for (int i = Coefficients.Length - 1; i >= 0; i--)
            {
                double c = Coefficients[i];
                if (c == 0.0)
                {
                    continue;
                }

                if (first)
                {
                    sb.Append(c.ToString(format, formatProvider));
                    if (i > 0)
                    {
                        sb.Append(VariableName);
                    }

                    if (i > 1)
                    {
                        sb.Append("^");
                        sb.Append(i);
                    }

                    first = false;
                }
                else
                {
                    if (c < 0.0)
                    {
                        sb.Append(" - ");
                        sb.Append((-c).ToString(format, formatProvider));
                    }
                    else
                    {
                        sb.Append(" + ");
                        sb.Append(c.ToString(format, formatProvider));
                    }

                    if (i > 0)
                    {
                        sb.Append(VariableName);
                    }

                    if (i > 1)
                    {
                        sb.Append("^");
                        sb.Append(i);
                    }
                }
            }

            return sb.ToString();
        }

        #endregion



        #region Equality

#pragma warning disable CS8767 // Nullability of reference types in type of parameter doesn't match implicitly implemented member (possibly because of nullability attributes).
        public bool Equals(Polynomial other)
#pragma warning restore CS8767 // Nullability of reference types in type of parameter doesn't match implicitly implemented member (possibly because of nullability attributes).
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;

            int n = Degree;
            if (n != other.Degree)
            {
                return false;
            }

            for (var i = 0; i <= n; i++)
            {
                if (!Coefficients[i].Equals(other.Coefficients[i]))
                {
                    return false;
                }
            }

            return true;
        }

#pragma warning disable CS8765 // Nullability of type of parameter doesn't match overridden member (possibly because of nullability attributes).
        public override bool Equals(object obj)
#pragma warning restore CS8765 // Nullability of type of parameter doesn't match overridden member (possibly because of nullability attributes).
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != typeof(Polynomial)) return false;
            return Equals((Polynomial)obj);
        }

        public override int GetHashCode()
        {
            var hashNum = Math.Min(Degree + 1, 25);
            int hash = 17;
            unchecked
            {
                for (var i = 0; i < hashNum; i++)
                {
                    hash = hash * 31 + Coefficients[i].GetHashCode();
                }
            }

            return hash;
        }

        #endregion

        #region Clone

        public Polynomial Clone()
        {
            int degree = EvaluateDegree(Coefficients);
            var coefficients = new double[degree + 1];
            Array.Copy(Coefficients, coefficients, coefficients.Length);
            return new Polynomial(coefficients);
        }

        object ICloneable.Clone()
        {
            return Clone();
        }

        #endregion
    }
}