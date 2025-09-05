using System;
using System.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using SpecialFunctions;

namespace SpecialFunctions.Tests
{
    [TestClass]
    public class SpecialFunctionsTests
    {
        private const double Tolerance = 1e-10;
        private const double RelativeTolerance = 1e-12;

        #region Helper Methods

        /// <summary>
        /// Assert that two doubles are approximately equal
        /// </summary>
        private void AssertApproximatelyEqual(double expected, double actual, double tolerance = Tolerance)
        {
            if (double.IsNaN(expected) && double.IsNaN(actual))
                return;
            
            if (double.IsInfinity(expected) && double.IsInfinity(actual))
            {
                Assert.AreEqual(Math.Sign(expected), Math.Sign(actual), "Infinity signs should match");
                return;
            }

            double diff = Math.Abs(expected - actual);
            double relativeDiff = Math.Abs(expected) > 1e-10 ? diff / Math.Abs(expected) : diff;
            
            Assert.IsTrue(diff < tolerance || relativeDiff < RelativeTolerance, 
                $"Expected: {expected}, Actual: {actual}, Difference: {diff}, Relative: {relativeDiff}");
        }

        /// <summary>
        /// Assert that two complex numbers are approximately equal
        /// </summary>
        private void AssertApproximatelyEqual(Complex expected, Complex actual, double tolerance = Tolerance)
        {
            AssertApproximatelyEqual(expected.Real, actual.Real, tolerance);
            AssertApproximatelyEqual(expected.Imaginary, actual.Imaginary, tolerance);
        }

        #endregion

        #region Basic Mathematical Functions Tests

        [TestMethod]
        public void Factorial_ValidInputs_ReturnsCorrectValues()
        {
            // Test known factorial values
            Assert.AreEqual(1, SF.factorial(0));
            Assert.AreEqual(1, SF.factorial(1));
            Assert.AreEqual(2, SF.factorial(2));
            Assert.AreEqual(6, SF.factorial(3));
            Assert.AreEqual(24, SF.factorial(4));
            Assert.AreEqual(120, SF.factorial(5));
            Assert.AreEqual(720, SF.factorial(6));
            Assert.AreEqual(5040, SF.factorial(7));
        }

        #endregion

        #region Airy Functions Tests

        [TestMethod]
        public void Airy_RealInputs_ReturnsExpectedValues()
        {
            // Test Airy function Ai(0) = 3^(-2/3)/Γ(2/3) ≈ 0.355028053887817
            double ai0 = SF.airy(0.0);
            AssertApproximatelyEqual(0.355028053887817, ai0, 1e-10);

            // Test Airy function at x = 1
            double ai1 = SF.airy(1.0);
            AssertApproximatelyEqual(0.135292416312881, ai1, 1e-10);

            // Test negative values
            double aiNeg1 = SF.airy(-1.0);
            AssertApproximatelyEqual(0.535560883292764, aiNeg1, 1e-10);
        }

        [TestMethod]
        public void Airy_ComplexInputs_ReturnsValidResults()
        {
            Complex z = new Complex(1.0, 0.5);
            Complex result = SF.airy(z);
            
            // Should return a finite complex number
            Assert.IsFalse(double.IsNaN(result.Real));
            Assert.IsFalse(double.IsNaN(result.Imaginary));
            Assert.IsFalse(double.IsInfinity(result.Real));
            Assert.IsFalse(double.IsInfinity(result.Imaginary));
        }

        [TestMethod]
        public void Airy_DerivativeTypes_ReturnsValidResults()
        {
            double x = 1.0;
            
            // Test all derivative types
            double ai = SF.airy(0, x);    // Ai(x)
            double aiPrime = SF.airy(1, x);   // Ai'(x)
            double bi = SF.airy(2, x);    // Bi(x)
            double biPrime = SF.airy(3, x);   // Bi'(x)
            
            // All should be finite
            Assert.IsFalse(double.IsNaN(ai));
            Assert.IsFalse(double.IsNaN(aiPrime));
            Assert.IsFalse(double.IsNaN(bi));
            Assert.IsFalse(double.IsNaN(biPrime));
            
            // Invalid type should return NaN
            Assert.IsTrue(double.IsNaN(SF.airy(4, x)));
        }

        #endregion

        #region Bessel Functions Tests

        [TestMethod]
        public void BesselJ_KnownValues_ReturnsCorrectResults()
        {
            // J₀(0) = 1
            AssertApproximatelyEqual(1.0, SF.besselj(0, 0.0));
            
            // J₁(0) = 0
            AssertApproximatelyEqual(0.0, SF.besselj(1, 0.0));
            
            // J₀(1) ≈ 0.7651976865579666
            AssertApproximatelyEqual(0.7651976865579666, SF.besselj(0, 1.0), 1e-10);
            
            // J₁(1) ≈ 0.4400505857449335
            AssertApproximatelyEqual(0.4400505857449335, SF.besselj(1, 1.0), 1e-10);
        }

        [TestMethod]
        public void BesselY_KnownValues_ReturnsCorrectResults()
        {
            // Y₀(1) ≈ 0.0882569642156769
            AssertApproximatelyEqual(0.0882569642156769, SF.bessely(0, 1.0), 1e-10);
            
            // Y₁(1) ≈ -0.7812128213002887
            AssertApproximatelyEqual(-0.7812128213002887, SF.bessely(1, 1.0), 1e-10);
        }

        [TestMethod]
        public void BesselI_KnownValues_ReturnsCorrectResults()
        {
            // I₀(0) = 1
            AssertApproximatelyEqual(1.0, SF.besseli(0, 0.0));
            
            // I₁(0) = 0
            AssertApproximatelyEqual(0.0, SF.besseli(1, 0.0));
            
            // I₀(1) ≈ 1.2660658777520083
            AssertApproximatelyEqual(1.2660658777520083, SF.besseli(0, 1.0), 1e-10);
        }

        [TestMethod]
        public void BesselK_KnownValues_ReturnsCorrectResults()
        {
            // K₀(1) ≈ 0.4210244382407083
            AssertApproximatelyEqual(0.4210244382407083, SF.besselk(0, 1.0), 1e-10);
            
            // K₁(1) ≈ 0.6019072301972346
            AssertApproximatelyEqual(0.6019072301972346, SF.besselk(1, 1.0), 1e-10);
        }

        [TestMethod]
        public void BesselH_ComplexInputs_ReturnsValidResults()
        {
            Complex z = new Complex(1.0, 0.5);
            
            Complex h1 = SF.besselh(0, z);
            Complex h2 = SF.besselh(0, 2, z);
            
            // Should return finite complex numbers
            Assert.IsFalse(double.IsNaN(h1.Real));
            Assert.IsFalse(double.IsNaN(h1.Imaginary));
            Assert.IsFalse(double.IsNaN(h2.Real));
            Assert.IsFalse(double.IsNaN(h2.Imaginary));
            
            // Invalid kind should return NaN
            Complex invalid = SF.besselh(0, 3, z);
            Assert.IsTrue(double.IsNaN(invalid.Real));
        }

        #endregion

        #region Spherical Bessel Functions Tests

        [TestMethod]
        public void SphericalBesselJ_KnownValues_ReturnsCorrectResults()
        {
            // j₀(0) = 1
            AssertApproximatelyEqual(1.0, SF.sbesselj(0, 0.0));
            
            // j₁(0) = 0
            AssertApproximatelyEqual(0.0, SF.sbesselj(1, 0.0));
            
            // j₀(π) = 0 (approximately)
            AssertApproximatelyEqual(0.0, SF.sbesselj(0, Math.PI), 1e-10);
        }

        [TestMethod]
        public void SphericalBesselJ_NegativeOrder_ReturnsNaN()
        {
            Assert.IsTrue(double.IsNaN(SF.sbesselj(-1, 1.0)));
        }

        [TestMethod]
        public void SphericalBesselY_KnownValues_ReturnsCorrectResults()
        {
            // y₀(π/2) should be finite
            double result = SF.sbessely(0, Math.PI / 2);
            Assert.IsFalse(double.IsNaN(result));
            Assert.IsFalse(double.IsInfinity(result));
        }

        [TestMethod]
        public void SphericalBesselY_AtZero_ReturnsNegativeInfinity()
        {
            Assert.IsTrue(double.IsNegativeInfinity(SF.sbessely(0, 0.0)));
        }

        #endregion

        #region Beta Functions Tests

        [TestMethod]
        public void Beta_KnownValues_ReturnsCorrectResults()
        {
            // B(1,1) = 1
            AssertApproximatelyEqual(1.0, SF.beta(1, 1));
            
            // B(2,1) = 1/2
            AssertApproximatelyEqual(0.5, SF.beta(2, 1));
            
            // B(1,2) = 1/2
            AssertApproximatelyEqual(0.5, SF.beta(1, 2));
            
            // B(2,2) = 1/6
            AssertApproximatelyEqual(1.0/6.0, SF.beta(2, 2), 1e-10);
        }

        [TestMethod]
        public void BetaLn_KnownValues_ReturnsCorrectResults()
        {
            // ln(B(1,1)) = ln(1) = 0
            AssertApproximatelyEqual(0.0, SF.betaln(1, 1));
            
            // ln(B(2,1)) = ln(1/2) = -ln(2)
            AssertApproximatelyEqual(-Math.Log(2), SF.betaln(2, 1));
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Beta_NegativeInput_ThrowsException()
        {
            SF.beta(-1, 1);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Beta_ZeroInput_ThrowsException()
        {
            SF.beta(0, 1);
        }

        #endregion

        #region Error Functions Tests

        [TestMethod]
        public void Erf_KnownValues_ReturnsCorrectResults()
        {
            // erf(0) = 0
            AssertApproximatelyEqual(0.0, SF.erf(0.0));
            
            // erf(∞) = 1
            Assert.AreEqual(1.0, SF.erf(double.PositiveInfinity));
            
            // erf(-∞) = -1
            Assert.AreEqual(-1.0, SF.erf(double.NegativeInfinity));
            
            // erf(1) ≈ 0.8427007929497149
            AssertApproximatelyEqual(0.8427007929497149, SF.erf(1.0), 1e-10);
            
            // erf is odd function: erf(-x) = -erf(x)
            AssertApproximatelyEqual(-SF.erf(1.0), SF.erf(-1.0));
        }

        [TestMethod]
        public void Erfc_KnownValues_ReturnsCorrectResults()
        {
            // erfc(0) = 1
            AssertApproximatelyEqual(1.0, SF.erfc(0.0));
            
            // erfc(∞) = 0
            Assert.AreEqual(0.0, SF.erfc(double.PositiveInfinity));
            
            // erfc(-∞) = 2
            Assert.AreEqual(2.0, SF.erfc(double.NegativeInfinity));
            
            // erfc(x) + erf(x) = 1
            double x = 1.5;
            AssertApproximatelyEqual(1.0, SF.erf(x) + SF.erfc(x));
        }

        [TestMethod]
        public void ErfInv_KnownValues_ReturnsCorrectResults()
        {
            // erfinv(0) = 0
            AssertApproximatelyEqual(0.0, SF.erfinv(0.0));
            
            // erfinv(erf(x)) = x for valid x
            double x = 0.5;
            AssertApproximatelyEqual(x, SF.erfinv(SF.erf(x)), 1e-10);
            
            // erfinv(1) = ∞
            Assert.AreEqual(double.PositiveInfinity, SF.erfinv(1.0));
            
            // erfinv(-1) = -∞
            Assert.AreEqual(double.NegativeInfinity, SF.erfinv(-1.0));
        }

        [TestMethod]
        public void ErfcInv_KnownValues_ReturnsCorrectResults()
        {
            // erfcinv(1) = 0
            AssertApproximatelyEqual(0.0, SF.erfcinv(1.0));
            
            // erfcinv(erfc(x)) = x for valid x
            double x = 0.5;
            AssertApproximatelyEqual(x, SF.erfcinv(SF.erfc(x)), 1e-10);
            
            // erfcinv(0) = ∞
            Assert.AreEqual(double.PositiveInfinity, SF.erfcinv(0.0));
            
            // erfcinv(2) = -∞
            Assert.AreEqual(double.NegativeInfinity, SF.erfcinv(2.0));
        }

        #endregion

        #region Legendre Functions Tests

        [TestMethod]
        public void Legendre0_KnownValues_ReturnsCorrectResults()
        {
            // P₀(x) = 1
            AssertApproximatelyEqual(1.0, SF.legendre0(0, 0.5));
            AssertApproximatelyEqual(1.0, SF.legendre0(0, -0.5));
            
            // P₁(x) = x
            AssertApproximatelyEqual(0.5, SF.legendre0(1, 0.5));
            AssertApproximatelyEqual(-0.5, SF.legendre0(1, -0.5));
            
            // P₂(x) = (3x² - 1)/2
            double x = 0.5;
            double expected = (3 * x * x - 1) / 2;
            AssertApproximatelyEqual(expected, SF.legendre0(2, x));
        }

        [TestMethod]
        public void Legendre_ArrayOutput_ReturnsCorrectLength()
        {
            int n = 5;
            double[] result = SF.legendre(n, 0.5);
            
            Assert.AreEqual(n + 1, result.Length);
            
            // First element should match P₀
            AssertApproximatelyEqual(SF.legendre0(n, 0.5), result[0]);
        }

        [TestMethod]
        public void Legendre_NormalizedVersion_ReturnsValidResults()
        {
            int n = 3;
            double x = 0.5;
            
            double[] normal = SF.legendre(n, x);
            double[] normalized = SF.legendre(n, x, "norm");
            double[] schmidt = SF.legendre(n, x, "sch");
            
            Assert.AreEqual(n + 1, normal.Length);
            Assert.AreEqual(n + 1, normalized.Length);
            Assert.AreEqual(n + 1, schmidt.Length);
            
            // Schmidt normalization should have P₀ = 0
            AssertApproximatelyEqual(0.0, schmidt[0]);
        }

        #endregion

        #region Gamma Functions Tests

        [TestMethod]
        public void Gamma_KnownValues_ReturnsCorrectResults()
        {
            // Γ(1) = 1
            AssertApproximatelyEqual(1.0, SF.gamma(1.0));
            
            // Γ(2) = 1
            AssertApproximatelyEqual(1.0, SF.gamma(2.0));
            
            // Γ(3) = 2
            AssertApproximatelyEqual(2.0, SF.gamma(3.0));
            
            // Γ(4) = 6
            AssertApproximatelyEqual(6.0, SF.gamma(4.0));
            
            // Γ(0.5) = √π
            AssertApproximatelyEqual(Math.Sqrt(Math.PI), SF.gamma(0.5), 1e-10);
        }

        [TestMethod]
        public void GammaLn_KnownValues_ReturnsCorrectResults()
        {
            // ln(Γ(1)) = ln(1) = 0
            AssertApproximatelyEqual(0.0, SF.gammaln(1.0));
            
            // ln(Γ(2)) = ln(1) = 0
            AssertApproximatelyEqual(0.0, SF.gammaln(2.0));
            
            // ln(Γ(3)) = ln(2)
            AssertApproximatelyEqual(Math.Log(2), SF.gammaln(3.0));
            
            // For large x, should be finite
            double largeLn = SF.gammaln(100.0);
            Assert.IsFalse(double.IsNaN(largeLn));
            Assert.IsFalse(double.IsInfinity(largeLn));
        }

        [TestMethod]
        public void GammaInc_ValidInputs_ReturnsFiniteResults()
        {
            double result1 = SF.gammainc(1.0, 2.0);
            double result2 = SF.gammainc(1.0, 2.0, "lower");
            double result3 = SF.gammainc(1.0, 2.0, "upper");
            
            Assert.IsFalse(double.IsNaN(result1));
            Assert.IsFalse(double.IsNaN(result2));
            Assert.IsFalse(double.IsNaN(result3));
            
            // Invalid type should return NaN
            Assert.IsTrue(double.IsNaN(SF.gammainc(1.0, 2.0, "invalid")));
        }

        [TestMethod]
        public void Psi_KnownValues_ReturnsCorrectResults()
        {
            // ψ(1) = -γ (Euler-Mascheroni constant ≈ -0.5772156649015329)
            AssertApproximatelyEqual(-0.5772156649015329, SF.psi(1.0), 1e-10);
            
            // ψ(2) = 1 - γ
            AssertApproximatelyEqual(1.0 - 0.5772156649015329, SF.psi(2.0), 1e-10);
        }

        [TestMethod]
        public void Psi_SpecialCases_HandlesCorrectly()
        {
            // ψ(0) should be -∞
            Assert.IsTrue(double.IsNegativeInfinity(SF.psi(0.0)));
            
            // ψ(-1) should be -∞  
            Assert.IsTrue(double.IsNegativeInfinity(SF.psi(-1.0)));
            
            // ψ(NaN) should be NaN
            Assert.IsTrue(double.IsNaN(SF.psi(double.NaN)));
        }

        #endregion

        #region Constants Tests

        [TestMethod]
        public void Constants_HaveExpectedValues()
        {
            // Test mathematical constants
            AssertApproximatelyEqual(Math.PI, SF.PI, 1e-10);
            AssertApproximatelyEqual(Math.PI / 180.0, SF.DEGRAD, 1e-15);
            
            // Test that geodetic constants are reasonable
            Assert.IsTrue(SF.a > 6e6 && SF.a < 7e6); // Earth radius in reasonable range
            Assert.IsTrue(SF.e2 > 0 && SF.e2 < 1);   // Eccentricity squared
            Assert.IsTrue(SF.Omega > 0);              // Angular velocity
            Assert.IsTrue(SF.GM > 0);                 // Gravitational parameter
        }

        #endregion

        #region Performance Tests (Optional)

        [TestMethod]
        public void Performance_GammaFunction_HandlesLargeInputs()
        {
            var stopwatch = System.Diagnostics.Stopwatch.StartNew();
            
            for (int i = 1; i <= 100; i++)
            {
                double result = SF.gamma(i * 0.1);
                Assert.IsFalse(double.IsNaN(result));
            }
            
            stopwatch.Stop();
            
            // Should complete in reasonable time (less than 1 second)
            Assert.IsTrue(stopwatch.ElapsedMilliseconds < 1000);
        }

        #endregion
    }

    #region Integration Tests

    [TestClass]
    public class SpecialFunctionsIntegrationTests
    {
        [TestMethod]
        public void Integration_ErrorFunctionIdentity_IsValid()
        {
            // Test that erf(x) + erfc(x) = 1 for various x values
            double[] testValues = { -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0 };
            
            foreach (double x in testValues)
            {
                double sum = SF.erf(x) + SF.erfc(x);
                Assert.AreEqual(1.0, sum, 1e-14, $"erf({x}) + erfc({x}) should equal 1");
            }
        }

        [TestMethod]
        public void Integration_BesselRecurrenceRelation_IsValid()
        {
            // Test Bessel recurrence relation: J_{n-1}(x) + J_{n+1}(x) = (2n/x) * J_n(x)
            double x = 2.0;
            int n = 1;
            
            double jPrev = SF.besselj(n - 1, x);  // J_0(2)
            double jCurr = SF.besselj(n, x);      // J_1(2)  
            double jNext = SF.besselj(n + 1, x);  // J_2(2)
            
            double leftSide = jPrev + jNext;
            double rightSide = (2.0 * n / x) * jCurr;
            
            Assert.AreEqual(leftSide, rightSide, 1e-12, "Bessel recurrence relation should hold");
        }
    }
    #endregion
}