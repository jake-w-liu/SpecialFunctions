using System;
using System.IO.Compression;
using Complex = System.Numerics.Complex;

namespace SpecialFunctions
{
    public static partial class SF
    {
        public static double PI = 3.1415926535897932384626430;
        public static double DEGRAD = PI / 180.0;
        public static double a = 6378137.0;
        public static double ep2 = 0.00673949677540;
        public static double Omega = 0.7292115E-4;
        public static double GM = 0.39860050000E+15;
        public static double e2 = 0.006694380022900;
        public static double k = 0.0019318513530;

        public const double TwoSqrtEOverPi = 1.8603827342052657173362492472666631120594218414085755;

        public const double SqrtPiOver2 = 1.2533141373155002512078826424055226265034933703050d;

        const int GammaN = 10;
        const double GammaR = 10.900511;

        static readonly double[] GammaDk =
        {
            2.48574089138753565546e-5,
            1.05142378581721974210,
            -3.45687097222016235469,
            4.51227709466894823700,
            -2.98285225323576655721,
            1.05639711577126713077,
            -1.95428773191645869583e-1,
            1.70970543404441224307e-2,
            -5.71926117404305781283e-4,
            4.63399473359905636708e-6,
            -2.71994908488607703910e-9
        };

        public static int factorial(int n)
        {
            int f = 1;
            if (n == 0)
            {
                return f;
            }
            else
            {
                f = n * factorial(n - 1);
            }
            return f;
        }

        public static Complex airy(Complex z)
        {
            return Amos.Cairy(z);
        }

        public static double airy(double z)
        {
            return airy(new Complex(z, 0)).Real;
        }

        public static Complex airy(int k, Complex z)
        {
            if (k == 1)
            {
                return Amos.CairyPrime(z);
            }
            else if (k == 2)
            {
                return Amos.Cbiry(z);
            }
            else if (k == 3)
            {
                return Amos.CbiryPrime(z);
            }
            else if (k == 0)
            {
                return airy(z);
            }
            else
            {
                return double.NaN;
            }
        }

        public static double airy(int k, double z)
        {
            Complex zc = new Complex(z, 0);
            return airy(k, zc).Real;
        }

        public static Complex besselj(double n, Complex z)
        {
            return Amos.Cbesj(n, z);
        }

        public static double besselj(double n, double z)
        {
            return Amos.Cbesj(n, z);
        }

        public static Complex bessely(double n, Complex z)
        {
            return Amos.Cbesy(n, z);
        }

        public static double bessely(double n, double z)
        {
            return Amos.Cbesy(n, z);
        }

        public static Complex besseli(double n, Complex z)
        {
            return Amos.Cbesi(n, z);
        }

        public static double besseli(double n, double z)
        {
            return besseli(n, new Complex(z, 0)).Real;
        }

        public static Complex besselk(double n, Complex z)
        {
            return Amos.Cbesk(n, z);
        }

        public static double besselk(double n, double z)
        {
            return Amos.Cbesk(n, z);
        }

        public static Complex besselh(double n, Complex z)
        {
            return Amos.Cbesh1(n, z);
        }

        public static Complex besselh(double n, int k, Complex z)
        {
            if (k == 2)
                return Amos.Cbesh2(n, z);
            else if (k == 1)
                return Amos.Cbesh1(n, z);
            else
                return double.NaN;
        }

        public static Complex sbesselj(double n, Complex z)
        {
            if (double.IsNaN(n) || double.IsNaN(z.Real) || double.IsNaN(z.Imaginary))
            {
                return new Complex(double.NaN, double.NaN);
            }

            if (double.IsInfinity(z.Real))
            {
                return (z.Imaginary == 0) ? Complex.Zero : new Complex(double.PositiveInfinity, double.PositiveInfinity);
            }

            if (z.Real == 0 && z.Imaginary == 0)
            {
                return (n == 0) ? 1 : 0;
            }

            return SqrtPiOver2 * besselj(n + 0.5, z) / Complex.Sqrt(z);
        }

        public static double sbesselj(double n, double z)
        {
            if (double.IsNaN(n) || double.IsNaN(z))
            {
                return double.NaN;
            }

            if (n < 0)
            {
                return double.NaN;
            }

            if (double.IsInfinity(z))
            {
                return 0;
            }

            if (z == 0)
            {
                return (n == 0) ? 1 : 0;
            }

            return SqrtPiOver2 * besselj(n + 0.5, z) / Math.Sqrt(z);
        }

        public static Complex sbessely(double n, Complex z)
        {
            if (double.IsNaN(n) || double.IsNaN(z.Real) || double.IsNaN(z.Imaginary))
            {
                return new Complex(double.NaN, double.NaN);
            }

            if (double.IsInfinity(z.Real))
            {
                return (z.Imaginary == 0) ? Complex.Zero : new Complex(double.PositiveInfinity, double.PositiveInfinity);
            }

            if (z.Real == 0 && z.Imaginary == 0)
            {
                return new Complex(double.NaN, double.NaN);
            }

            return SqrtPiOver2 * bessely(n + 0.5, z) / Complex.Sqrt(z);
        }

        public static double sbessely(double n, double z)
        {
            if (double.IsNaN(n) || double.IsNaN(z))
            {
                return double.NaN;
            }

            if (n < 0)
            {
                return double.NaN;
            }

            if (double.IsInfinity(z))
            {
                return 0;
            }

            if (z == 0)
            {
                return double.NegativeInfinity;
            }

            return SqrtPiOver2 * bessely(n + 0.5, z) / Math.Sqrt(z);
        }

        public static Complex sbesselh(double n, Complex z)
        {
            return sbesselj(n, z) + Complex.ImaginaryOne * sbessely(n, z);
        }

        public static Complex sbesselh(double n, int k, Complex z)
        {
            if (k == 1)
                return sbesselj(n, z) + Complex.ImaginaryOne * sbessely(n, z);
            else if (k == 2)
                return sbesselj(n, z) - Complex.ImaginaryOne * sbessely(n, z);
            else
                return Complex.NaN;
        }

        public static double betaln(double z, double w)
        {
            if (z <= 0.0)
            {
                throw new ArgumentException("Value must be positive.", nameof(z));
            }

            if (w <= 0.0)
            {
                throw new ArgumentException("Value must be positive.", nameof(w));
            }

            return gammaln(z) + gammaln(w) - gammaln(z + w);
        }

        public static double beta(double z, double w)
        {
            return Math.Exp(betaln(z, w));
        }

        public static double erf(double x)
        {
            if (x == 0)
            {
                return 0;
            }

            if (double.IsPositiveInfinity(x))
            {
                return 1;
            }

            if (double.IsNegativeInfinity(x))
            {
                return -1;
            }

            if (double.IsNaN(x))
            {
                return double.NaN;
            }

            return Amos.ErfImp(x, false);
        }

        public static double erfc(double x)
        {
            if (x == 0)
            {
                return 1;
            }

            if (double.IsPositiveInfinity(x))
            {
                return 0;
            }

            if (double.IsNegativeInfinity(x))
            {
                return 2;
            }

            if (double.IsNaN(x))
            {
                return double.NaN;
            }

            return Amos.ErfImp(x, true);
        }

        public static double erfinv(double z)
        {
            if (z == 0.0)
            {
                return 0.0;
            }

            if (z >= 1.0)
            {
                return double.PositiveInfinity;
            }

            if (z <= -1.0)
            {
                return double.NegativeInfinity;
            }

            double p, q, s;
            if (z < 0)
            {
                p = -z;
                q = 1 - p;
                s = -1;
            }
            else
            {
                p = z;
                q = 1 - z;
                s = 1;
            }

            return Amos.ErfInvImpl(p, q, s);
        }

        public static double erfcinv(double z)
        {
            if (z <= 0.0)
            {
                return double.PositiveInfinity;
            }

            if (z >= 2.0)
            {
                return double.NegativeInfinity;
            }

            double p, q, s;
            if (z > 1)
            {
                q = 2 - z;
                p = 1 - q;
                s = -1;
            }
            else
            {
                p = 1 - z;
                q = z;
                s = 1;
            }

            return Amos.ErfInvImpl(p, q, s);
        }

        public static double legendre0(int n, double x)
        {
            double LFS0 = 1.0;
            double LFS1 = x;
            double Lf = 0.0;
            if (n == 0)
            {
                Lf = LFS0;
            }
            else if (n == 1)
            {
                Lf = LFS1;
            }
            else if (n >= 1)
            {
                for (int i = 2; i <= n; i++)
                {
                    Lf = (-(i - 1) * LFS0 + (2 * i - 1) * x * LFS1) / i;
                    LFS0 = LFS1;
                    LFS1 = Lf;
                }
            }
            return Lf;
        }

        public static double[] legendre(int n, double x)
        {
            double[] Pnm = new double[n + 1];
            double[] ALFS = new double[n + 1];
            double[] BLFS = new double[n + 1];

            double s = Math.Sqrt(1.0 - x * x);

            switch (n)
            {
                case 0:
                    Pnm[0] = 1.0;
                    break;
                case 1:
                    Pnm[0] = x;
                    Pnm[1] = s;
                    break;
                case 2:
                    Pnm[0] = (3.0 * x * x - 1.0) / 2.0;
                    Pnm[1] = 3.0 * s * x;
                    Pnm[2] = 3.0 * (1.0 - x * x);
                    break;
                case 3:
                    Pnm[0] = x * (5.0 * x * x - 3.0) / 2.0;
                    Pnm[1] = 3.0 * (5.0 * x * x - 1.0) * s / 2.0;
                    Pnm[2] = 15.0 * x * (1.0 - x * x);
                    Pnm[3] = 15.0 * s * s * s;
                    break;
                case 4:
                    Pnm[0] = (35.0 * x * x * x * x - 30.0 * x * x + 3.0) / 8.0;
                    Pnm[1] = 5.0 * (7.0 * x * x - 3.0) * x * s / 2.0;
                    Pnm[2] = 15.0 * (7.0 * x * x - 1.0) * (1.0 - x * x) / 2.0;
                    Pnm[3] = 105.0 * s * s * s * x;
                    Pnm[4] = 105.0 * s * s * s * s;
                    break;
                case 5:
                    Pnm[0] = x * (63.0 * x * x * x * x - 70.0 * x * x + 15.0) / 8.0;
                    Pnm[1] = 15.0 * s * (21.0 * x * x * x * x - 14.0 * x * x + 1.0) / 8.0;
                    Pnm[2] = 105.0 * x * (3.0 * x * x - 1.0) * (1.0 - x * x) / 2.0;
                    Pnm[3] = 105.0 * s * s * s * (9.0 * x * x - 1.0) / 2.0;
                    Pnm[4] = 945.0 * s * s * s * s * x;
                    Pnm[5] = 945.0 * s * s * s * s * s;
                    break;
                case 6:
                    Pnm[0] = (x * x * (x * x * (231.0 * x * x - 315.0) + 105.0) - 5.0) / 16.0;
                    Pnm[1] = 21.0 * x * (x * x * (33.0 * x * x - 30.0) + 5.0) * s / 8.0;
                    Pnm[2] = 105.0 * s * s * (x * x * (33.0 * x * x - 18.0) + 1.0) / 8.0;
                    Pnm[3] = 315.0 * (11.0 * x * x - 3.0) * x * s * s * s / 2.0;
                    Pnm[4] = 945.0 * s * s * s * s * (11.0 * x * x - 1.0) / 2.0;
                    Pnm[5] = 10395.0 * x * s * s * s * s * s;
                    Pnm[6] = 10395.0 * s * s * s * s * s * s;
                    break;
                case 7:
                    Pnm[0] = x * (x * x * (429.0 * Math.Pow(x, 4) - 693.0 * x * x + 315.0) - 35.0) / 16.0;
                    Pnm[1] = 7.0 * s * (x * x * (429.0 * Math.Pow(x, 4) - 495.0 * x * x + 135.0) - 5.0) / 16.0;
                    Pnm[2] = 63.0 * x * s * s * (x * x * (143.0 * x * x - 110.0) + 15.0) / 8.0;
                    Pnm[3] = 315.0 * s * s * s * (x * x * (143.0 * x * x - 66.0) + 3.0) / 8.0;
                    Pnm[4] = 3465.0 * x * s * s * s * s * (13.0 * x * x - 3.0) / 2.0;
                    Pnm[5] = 10395.0 * Math.Pow(s, 5) * (13.0 * x * x - 1.0) / 2.0;
                    Pnm[6] = 135135.0 * x * Math.Pow(s, 6);
                    Pnm[7] = 135135.0 * Math.Pow(s, 7);
                    break;
                default:
                    ALFS[0] = (x * x * (x * x * (x * x * (6435.0 * x * x - 12012.0) + 6930.0) - 1260.0) + 35) / 128.0;
                    ALFS[1] = 9.0 * x * s * (x * x * (x * x * (715.0 * x * x - 1001.0) + 385.0) - 35.0) / 16.0;
                    ALFS[2] = 315.0 * s * s * (x * x * (x * x * (143.0 * x * x - 143.0) + 33.0) - 1.0) / 16.0;
                    ALFS[3] = 3465.0 * x * s * s * s * (x * x * (39.0 * x * x - 26.0) + 3.0) / 8.0;
                    ALFS[4] = 10395.0 * s * s * s * s * (65.0 * Math.Pow(x, 4) - 26.0 * x * x + 1.0) / 8.0;
                    ALFS[5] = 135135.0 * x * Math.Pow(s, 5) * (5.0 * x * x - 1.0) / 2.0;
                    ALFS[6] = 135135.0 * Math.Pow(s, 6) * (15.0 * x * x - 1.0) / 2.0;
                    ALFS[7] = 2027025.0 * x * s * s * s * s * s * s * s;
                    ALFS[8] = 2027025.0 * s * s * s * s * s * s * s * s;

                    if (n == 8)
                    {
                        Pnm = ALFS;
                    }
                    else
                    {
                        BLFS[0] = x * (x * x * (429.0 * Math.Pow(x, 4) - 693.0 * x * x + 315.0) - 35.0) / 16.0;
                        BLFS[1] = 7.0 * s * (x * x * (429.0 * Math.Pow(x, 4) - 495.0 * x * x + 135.0) - 5.0) / 16.0;
                        BLFS[2] = 63.0 * x * s * s * (x * x * (143.0 * x * x - 110.0) + 15.0) / 8.0;
                        BLFS[3] = 315.0 * s * s * s * (x * x * (143.0 * x * x - 66.0) + 3.0) / 8.0;
                        BLFS[4] = 3465.0 * x * s * s * s * s * (13.0 * x * x - 3.0) / 2.0;
                        BLFS[5] = 10395.0 * Math.Pow(s, 5) * (13.0 * x * x - 1.0) / 2.0;
                        BLFS[6] = 135135.0 * x * Math.Pow(s, 6);
                        BLFS[7] = 135135.0 * Math.Pow(s, 7);

                        int l = 8;
                        while (l < n)
                        {
                            Pnm[0] = legendre0(l + 1, x);
                            int k = 1;
                            while (k <= l + 1)
                            {
                                if (k <= l - 1)
                                {
                                    Pnm[k] = ((2 * l + 1) * x * ALFS[k] - (l + k) * BLFS[k]) / (l - k + 1);
                                }

                                else if (k > l - 1)
                                {
                                    Pnm[k] = -((l - k + 2) * x * Pnm[k - 1] - (l + k) * ALFS[k - 1]) / s;
                                }

                                k = k + 1;
                            }
                            Array.Copy(ALFS, BLFS, ALFS.Length);
                            Array.Copy(Pnm, ALFS, Pnm.Length);
                            l = l + 1;
                        }
                    }
                    break;
            }
            for (int I = 0; I < Pnm.Length; I++)
            {
                if (I % 2 == 1)
                    Pnm[I] = -Pnm[I];
            }
            return Pnm;
        }

        public static double[] legendre(int n, double x, string s)
        {
            double[] Pnm = new double[n + 1];
            Pnm = legendre(n, x);
            if (s == "nrom")
            {
                double logFactor;
                double normalization;
                double sign;
                for (int m = 0; m < Pnm.Length; m++)
                {
                    logFactor = Math.Log(n + 0.5)
                     + gammaln(n - m + 1)
                     - gammaln(n + m + 1);
                    normalization = Math.Exp(0.5 * logFactor);
                    sign = (m % 2 == 0) ? 1.0 : -1.0;
                    Pnm[m] = normalization * Pnm[m] * sign;
                    // Pnm[m] = Math.Sqrt((n + 0.5) * factorial(n - m) / factorial(n + m)) * Pnm[m] * Math.Pow(-1, m);
                }
            }
            else if (s == "sch")
            {
                double logFactor;
                double normalization;
                double sign;
                for (int m = 0; m < Pnm.Length; m++)
                {
                    logFactor = Math.Log(2)
                     + gammaln(n - m + 1)
                     - gammaln(n + m + 1);
                    normalization = Math.Exp(0.5 * logFactor);
                    sign = (m % 2 == 0) ? 1.0 : -1.0;
                    Pnm[m] = normalization * Pnm[m] * sign;
                }
                Pnm[0] = 0.0;
            }
            return Pnm;
        }

        public static double gamma(double z)
        {
            if (z < 0.5)
            {
                double s = GammaDk[0];
                for (int i = 1; i <= GammaN; i++)
                {
                    s += GammaDk[i] / (i - z);
                }

                return Math.PI / (Math.Sin(Math.PI * z)
                                * s
                                * TwoSqrtEOverPi
                                * Math.Pow((0.5 - z + GammaR) / Math.E, 0.5 - z));
            }
            else
            {
                double s = GammaDk[0];
                for (int i = 1; i <= GammaN; i++)
                {
                    s += GammaDk[i] / (z + i - 1.0);
                }

                return s * TwoSqrtEOverPi * Math.Pow((z - 0.5 + GammaR) / Math.E, z - 0.5);
            }
        }

        public static double gammainc(double x, double a)
        {
            return Amos.GammaLowerRegularized(a, x) * gamma(a);
        }

        public static double gammainc(double x, double a, string t)
        {
            if (t == "lower")
                return Amos.GammaLowerRegularized(a, x) * gamma(a);
            else if (t == "upper")
                return Amos.GammaLowerRegularized(a, x) * gamma(a);
            else
                return double.NaN;
        }

        public static double gammaln(double x)
        {
            // Lanczos approximation coefficients (g=7, n=9)
            double[] p = {
            0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7
        };

            if (x < 0.5)
            {
                // Reflection formula for x < 0.5
                return Math.Log(PI) - Math.Log(Math.Sin(PI * x)) - gammaln(1 - x);
            }

            x -= 1;
            double a = p[0];
            double t = x + 7.5; // g + 0.5
            for (int i = 1; i < p.Length; i++)
            {
                a += p[i] / (x + i);
            }

            return 0.5 * Math.Log(2 * PI) + (x + 0.5) * Math.Log(t) - t + Math.Log(a);
        }
        
        public static double psi(double x)
        {
            const double c = 12.0;
            const double d1 = -0.57721566490153286;
            const double d2 = 1.6449340668482264365;
            const double s = 1e-6;
            const double s3 = 1.0/12.0;
            const double s4 = 1.0/120.0;
            const double s5 = 1.0/252.0;
            const double s6 = 1.0/240.0;
            const double s7 = 1.0/132.0;

            if (double.IsNegativeInfinity(x) || double.IsNaN(x))
            {
                return double.NaN;
            }

            // Handle special cases.
            if (x <= 0 && Math.Floor(x) == x)
            {
                return double.NegativeInfinity;
            }

            // Use inversion formula for negative numbers.
            if (x < 0)
            {
                return psi(1.0 - x) + (Math.PI/Math.Tan(-Math.PI*x));
            }

            if (x <= s)
            {
                return d1 - (1/x) + (d2*x);
            }

            double result = 0;
            while (x < c)
            {
                result -= 1/x;
                x++;
            }

            if (x >= c)
            {
                var r = 1/x;
                result += Math.Log(x) - (0.5*r);
                r *= r;

                result -= r*(s3 - (r*(s4 - (r*(s5 - (r*(s6 - (r*s7))))))));
            }
            return result;
        }
    }
}