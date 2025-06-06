
# SpecialFunctions

A minimal C# library providing numerical implementations of common special functions, such as:

- Airy functions: Ai(z), Ai'(z), Bi(z), Bi'(z)
- Bessel functions: Jₙ(z), Yₙ(z), Iₙ(z), Kₙ(z), Hₙ^{(1)}(z), Hₙ^{(2)}(z)
- Spherical Bessel functions: jₙ(z), yₙ(z), hₙ^{(1)}(z), hₙ^{(2)}(z)
- Error functions: erf(x), erfc(x), erfinv(x), erfcinv(x)
- Beta and incomplete beta functions: beta(z, w), betaln(z, w)
- Legendre polynomials: Pₙ(x) and associated values
- Factorials and constants (π, Earth constants, etc.)

The package is built upon [Math.NET Numerics](https://github.com/mathnet/mathnet-numerics) and [AminKH's Legendre Polynomials repository](https://github.com/AminKH/Legendre-Polynomials). The API is designed to closely follow the structure and naming conventions of MATLAB, facilitating a smooth transition for users who are familiar with MATLAB’s syntax and functionality. This alignment allows for rapid prototyping, intuitive usage, and easier cross-platform adaptation of existing codebases.

## Usage

```csharp
using static SpecialFunctions.SF;
using System.Numerics;

double result = erf(1.0);
Complex bessel = besselj(1.5, new Complex(1.0, 0.5));

```

For detailed functionalities, please visit the [manual page](./docs/MANUAL.md).

## Notes

- Most functions are overloaded for both real (`double`) and complex (`System.Numerics.Complex`) inputs.
- Many implementations delegate to the `Amos` library for accurate numerical evaluation.
- Ensure input values are valid to avoid `NaN` results.

## License

MIT License or as applicable.
