# SmartBigRational: a big rational type optimized for addition, subtraction and multiplication

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
