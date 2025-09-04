---
title: 'SpecialFunctions: A C# package of special functions for scientific computing with MATLAB-compatible API'
tags:
  - C#
  - special functions
  - scientific computing
authors:
  - name: Jake W. Liu
    affiliation: 1
    corresponding: true
    orcid: 0000-0001-5458-7917
affiliations:
  - name: Department of Electronic Engineering, National Taipei University of Technology, Taiwan
    index: 1
date: 14 June 2025
bibliography: paper.bib
---

# Summary

`SpecialFunctions` is a C# library designed for scientific computing applications requiring a broad range of special functions, such as those found in mathematical physics and engineering [@abramowitz1965handbook; @arfken2011mathematical]. While Math.NET Numerics is the dominant open-source library for numerical computing in C#, its support for special functions remains incomplete. Furthermore, alternatives like ALGLIB are restricted by GPL licensing, making them unsuitable for commercial use.

This package addresses these gaps by offering a more complete set of special functions under the permissive MIT license, facilitating both academic and commercial adoption. The implementation leverages Math.NET Numerics where applicable [@mathnet] and incorporates extensions such as Legendre functions based on open-source contributions (e.g., AminKH's repository) [@lp]. The API design follows MATLAB conventions to support easier transition for users familiar with MATLAB, and includes functions like Bessel functions (including spherical variants), Airy functions, incomplete gamma functions, error functions, Legendre polynomials, and logarithmic derivatives.

The library is developed to support electromagnetic simulations and field transformations—specifically, cylindrical and spherical near-to-far-field (NTFF) calculations—in environments where C# is frequently used for hardware control. Full documentation and MATLAB-style references are provided for user guidance.

# Statement of need

The C# ecosystem lacks a comprehensive, permissively licensed library that offers a wide range of special functions required in scientific and engineering computations. Existing libraries such as Math.NET Numerics provide partial coverage but omit important function families such as associated Legendre polynomials. Other alternatives like ALGLIB include more functionality but are distributed under the GPL license, restricting their use in commercial applications.

This project addresses that gap by providing an open-source, MIT-licensed library that combines ease of use with broad coverage of special functions. It is designed to support both academic research and industrial development, particularly in domains where C# is widely adopted, such as real-time control systems and electromagnetic simulation software. The implementation adheres to MATLAB-style API conventions, facilitating a low-friction transition for users familiar with numerical computing in MATLAB.

By offering consistent naming conventions, support for complex-valued inputs, and compatibility with Math.NET Numerics, this library serves as a practical and extensible foundation for scientific applications—including, but not limited to, near-to-far-field (NTFF) transformations, wave propagation modeling, and geodetic calculations.

# Overview of features

The `SpecialFunctions` library provides a comprehensive set of special functions commonly used in scientific computing. The implementation is structured to be extensible, allowing for future inclusion of additional functions. Key features include:

- **Complete support for Bessel functions**:
  - Regular Bessel functions: `besselj`, `bessely`, `besseli`, `besselk`
  - Hankel functions: `besselh`
  - Spherical Bessel functions: `sbesselj`, `sbessely`, `sbesselh`

- **Airy functions**:
  - `airy(z)` and higher-order variants `airy(k, z)` for real and complex inputs

- **Gamma and related functions**:
  - Gamma function `gamma`
  - Log-Gamma function `gammaln`
  - Incomplete gamma function `gammainc` (lower and upper)
  - Digamma (Psi) function `psi`

- **Error functions**:
  - Error function `erf`
  - Complementary error function `erfc`
  - Inverse error functions: `erfinv`, `erfcinv`

- **Legendre polynomials**:
  - Legendre polynomials: `legendre0(n, x)`
  - Associated Legendre functions with normalization options: `legendre(n, x, mode)`
    - Supports `"nrom"` (normalized) and `"sch"` (Schmidt semi-normalized) modes

- **Beta function and logarithmic form**:
  - `beta(z, w)`, `betaln(z, w)`

Each function supports real and complex input variants where applicable. The library is implemented as a static class (`SF`) under the `SpecialFunctions` namespace.

# Acknowledgements

Legendre function implementations reference AminKH’s repository. The library also relies on selected components from Math.NET Numerics.

# References
