# Special Functions Library API Documentation

A C# library providing special mathematical functions with MATLAB-compatible naming conventions.

## Namespace
```csharp
SpecialFunctions.SF
```

## Functions

### Basic Mathematical Functions

#### `factorial`
```csharp
public static int factorial(int n)
```
Computes the factorial of a non-negative integer.

**Parameters:**
- `n` (int): Non-negative integer

**Returns:** 
- `int`: n! = n × (n-1) × ... × 1

**Example:**
```csharp
int result = SF.factorial(5); // Returns 120
```

---

### Airy Functions

#### `airy`
```csharp
public static Complex airy(Complex z)
public static double airy(double z)
public static Complex airy(int k, Complex z)
public static double airy(int k, double z)
```
Computes Airy functions and their derivatives.

**Parameters:**
- `z` (Complex/double): Input argument
- `k` (int): Function type selector
  - `0`: Ai(z) - Airy function of the first kind
  - `1`: Ai'(z) - Derivative of Airy function of the first kind
  - `2`: Bi(z) - Airy function of the second kind
  - `3`: Bi'(z) - Derivative of Airy function of the second kind

**Returns:**
- `Complex/double`: Value of the specified Airy function

**Example:**
```csharp
double ai = SF.airy(1.0);           // Ai(1)
double aiPrime = SF.airy(1, 1.0);   // Ai'(1)
Complex aiComplex = SF.airy(new Complex(1, 0.5));
```

---

### Bessel Functions

#### `besselj`
```csharp
public static Complex besselj(double n, Complex z)
public static double besselj(double n, double z)
```
Computes Bessel functions of the first kind Jₙ(z).

**Parameters:**
- `n` (double): Order of the Bessel function
- `z` (Complex/double): Input argument

**Returns:**
- `Complex/double`: Jₙ(z)

#### `bessely`
```csharp
public static Complex bessely(double n, Complex z)
public static double bessely(double n, double z)
```
Computes Bessel functions of the second kind Yₙ(z).

**Parameters:**
- `n` (double): Order of the Bessel function
- `z` (Complex/double): Input argument

**Returns:**
- `Complex/double`: Yₙ(z)

#### `besseli`
```csharp
public static Complex besseli(double n, Complex z)
public static double besseli(double n, double z)
```
Computes modified Bessel functions of the first kind Iₙ(z).

**Parameters:**
- `n` (double): Order of the Bessel function
- `z` (Complex/double): Input argument

**Returns:**
- `Complex/double`: Iₙ(z)

#### `besselk`
```csharp
public static Complex besselk(double n, Complex z)
public static double besselk(double n, double z)
```
Computes modified Bessel functions of the second kind Kₙ(z).

**Parameters:**
- `n` (double): Order of the Bessel function
- `z` (Complex/double): Input argument

**Returns:**
- `Complex/double`: Kₙ(z)

#### `besselh`
```csharp
public static Complex besselh(double n, Complex z)
public static Complex besselh(double n, int k, Complex z)
```
Computes Hankel functions (Bessel functions of the third kind).

**Parameters:**
- `n` (double): Order of the Bessel function
- `z` (Complex): Input argument
- `k` (int): Kind selector (1 or 2)

**Returns:**
- `Complex`: Hₙ⁽¹⁾(z) or Hₙ⁽²⁾(z)

**Example:**
```csharp
Complex h1 = SF.besselh(0.5, new Complex(2, 1));     // H₀.₅⁽¹⁾(2+i)
Complex h2 = SF.besselh(0.5, 2, new Complex(2, 1));  // H₀.₅⁽²⁾(2+i)
```

---

### Spherical Bessel Functions

#### `sbesselj`
```csharp
public static Complex sbesselj(double n, Complex z)
public static double sbesselj(double n, double z)
```
Computes spherical Bessel functions of the first kind jₙ(z).

**Parameters:**
- `n` (double): Order (must be non-negative for real inputs)
- `z` (Complex/double): Input argument

**Returns:**
- `Complex/double`: jₙ(z) = √(π/2z) × Jₙ₊₀.₅(z)

#### `sbessely`
```csharp
public static Complex sbessely(double n, Complex z)
public static double sbessely(double n, double z)
```
Computes spherical Bessel functions of the second kind yₙ(z).

**Parameters:**
- `n` (double): Order (must be non-negative for real inputs)
- `z` (Complex/double): Input argument

**Returns:**
- `Complex/double`: yₙ(z) = √(π/2z) × Yₙ₊₀.₅(z)

#### `sbesselh`
```csharp
public static Complex sbesselh(double n, Complex z)
public static Complex sbesselh(double n, int k, Complex z)
```
Computes spherical Hankel functions.

**Parameters:**
- `n` (double): Order
- `z` (Complex): Input argument
- `k` (int): Kind selector (1 or 2)

**Returns:**
- `Complex`: hₙ⁽¹⁾(z) or hₙ⁽²⁾(z)

---

### Beta Functions

#### `beta`
```csharp
public static double beta(double z, double w)
```
Computes the beta function B(z,w).

**Parameters:**
- `z` (double): First parameter (must be positive)
- `w` (double): Second parameter (must be positive)

**Returns:**
- `double`: B(z,w) = Γ(z)Γ(w)/Γ(z+w)

#### `betaln`
```csharp
public static double betaln(double z, double w)
```
Computes the natural logarithm of the beta function.

**Parameters:**
- `z` (double): First parameter (must be positive)
- `w` (double): Second parameter (must be positive)

**Returns:**
- `double`: ln(B(z,w))

**Example:**
```csharp
double b = SF.beta(2.5, 3.0);
double logB = SF.betaln(2.5, 3.0);
```

---

### Error Functions

#### `erf`
```csharp
public static double erf(double x)
```
Computes the error function erf(x).

**Parameters:**
- `x` (double): Input value

**Returns:**
- `double`: erf(x) = (2/√π) ∫₀ˣ e⁻ᵗ² dt

#### `erfc`
```csharp
public static double erfc(double x)
```
Computes the complementary error function erfc(x).

**Parameters:**
- `x` (double): Input value

**Returns:**
- `double`: erfc(x) = 1 - erf(x)

#### `erfinv`
```csharp
public static double erfinv(double z)
```
Computes the inverse error function.

**Parameters:**
- `z` (double): Input value (-1 < z < 1)

**Returns:**
- `double`: y such that erf(y) = z

#### `erfcinv`
```csharp
public static double erfcinv(double z)
```
Computes the inverse complementary error function.

**Parameters:**
- `z` (double): Input value (0 < z < 2)

**Returns:**
- `double`: y such that erfc(y) = z

**Example:**
```csharp
double erfVal = SF.erf(1.0);           // ≈ 0.8427
double erfcVal = SF.erfc(1.0);         // ≈ 0.1573
double invErf = SF.erfinv(0.5);        // ≈ 0.4769
```

---

### Legendre Functions

#### `legendre0`
```csharp
public static double legendre0(int n, double x)
```
Computes Legendre polynomial Pₙ(x).

**Parameters:**
- `n` (int): Degree of the polynomial
- `x` (double): Input value (-1 ≤ x ≤ 1)

**Returns:**
- `double`: Pₙ(x)

#### `legendre`
```csharp
public static double[] legendre(int n, double x)
public static double[] legendre(int n, double x, string s)
```
Computes associated Legendre functions Pₙᵐ(x) for m = 0, 1, ..., n.

**Parameters:**
- `n` (int): Maximum degree
- `x` (double): Input value (-1 ≤ x ≤ 1)
- `s` (string): Normalization type (optional)
  - `"nrom"`: Normalized Legendre functions
  - `"sch"`: Schmidt semi-normalized functions

**Returns:**
- `double[]`: Array of length (n+1) containing [P₀ⁿ(x), P₁ⁿ(x), ..., Pₙⁿ(x)]

**Example:**
```csharp
double p3 = SF.legendre0(3, 0.5);              // P₃(0.5)
double[] pnm = SF.legendre(3, 0.5);             // [P₀³, P₁³, P₂³, P₃³]
double[] pnmNorm = SF.legendre(3, 0.5, "nrom"); // Normalized functions
```

---

### Gamma Functions

#### `gamma`
```csharp
public static double gamma(double z)
```
Computes the gamma function Γ(z).

**Parameters:**
- `z` (double): Input value

**Returns:**
- `double`: Γ(z)

#### `gammaln`
```csharp
public static double gammaln(double x)
```
Computes the natural logarithm of the gamma function.

**Parameters:**
- `x` (double): Input value

**Returns:**
- `double`: ln(Γ(x))

#### `gammainc`
```csharp
public static double gammainc(double x, double a)
public static double gammainc(double x, double a, string t)
```
Computes the incomplete gamma function.

**Parameters:**
- `x` (double): Upper limit of integration
- `a` (double): Parameter
- `t` (string): Type ("lower" or "upper")

**Returns:**
- `double`: Incomplete gamma function value

#### `psi`
```csharp
public static double psi(double x)
```
Computes the digamma function ψ(x) (derivative of ln(Γ(x))).

**Parameters:**
- `x` (double): Input value

**Returns:**
- `double`: ψ(x) = Γ'(x)/Γ(x)

**Example:**
```csharp
double g = SF.gamma(2.5);              // Γ(2.5)
double logG = SF.gammaln(2.5);         // ln(Γ(2.5))
double incG = SF.gammainc(1.0, 2.0);   // Incomplete gamma
double psiVal = SF.psi(2.0);           // ψ(2)
```

---
