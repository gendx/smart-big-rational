# SmartBigRational: a big rational type optimized for addition, subtraction and multiplication

[![Crate](https://img.shields.io/crates/v/smart-big-rational.svg?logo=rust)](https://crates.io/crates/smart-big-rational)
[![Documentation](https://img.shields.io/docsrs/smart-big-rational/0.1.0?logo=rust)](https://docs.rs/smart-big-rational/0.1.0/)
[![Minimum Rust 1.87.0](https://img.shields.io/crates/msrv/smart-big-rational/0.1.0.svg?logo=rust&color=orange)](https://releases.rs/docs/1.87.0/)
[![Lines of Code](https://www.aschey.tech/tokei/github/gendx/smart-big-rational?category=code&branch=0.1.0)](https://github.com/gendx/smart-big-rational/tree/0.1.0)
[![Dependencies](https://deps.rs/crate/smart-big-rational/0.1.0/status.svg)](https://deps.rs/crate/smart-big-rational/0.1.0)
[![License](https://img.shields.io/crates/l/smart-big-rational/0.1.0.svg)](https://github.com/gendx/smart-big-rational/blob/0.1.0/LICENSE)
[![Codecov](https://codecov.io/gh/gendx/smart-big-rational/branch/0.1.0/graph/badge.svg)](https://codecov.io/gh/gendx/smart-big-rational/tree/0.1.0)
[![Build Status](https://github.com/gendx/smart-big-rational/actions/workflows/build.yml/badge.svg?branch=0.1.0)](https://github.com/gendx/smart-big-rational/actions/workflows/build.yml)
[![Test Status](https://github.com/gendx/smart-big-rational/actions/workflows/tests.yml/badge.svg?branch=0.1.0)](https://github.com/gendx/smart-big-rational/actions/workflows/tests.yml)

This Rust library provides the `SmartBigRational` type, similar to the `num`
crate's
[`BigRational`](https://docs.rs/num-rational/latest/num_rational/type.BigRational.html)
type, but optimized for addition, subtraction and multiplication operations.

Under the hood, the denominator is represented as a product of small primes,
which allows fast Greatest Common Divisor (GCD) operations, which are often the
bottleneck with the vanilla `BigRational` type. Additionally, `SmartBigRational`
aren't always represented as reduced fractions, which again avoids GCD
operations.

The underlying technique is described in the blog post
[_Optimization adventures: making a parallel Rust workload even faster with data-oriented design (and other tricks)_](https://gendx.dev/blog/2024/12/02/rust-data-oriented-design.html#optimizing-bigrationals-15x-to-20x-faster).

On the flip side, memory usage with non-reduced `SmartBigRational` fractions can
increase, notably when computing an expression with high multiplicative depth.
In these cases, vanilla `BigRational` may be more efficient.

In practice, `SmartBigRational` are well suited for expressions with long sums,
low multiplicative depth and/or small denominators.

As always: benchmark before choosing one or the other option.
