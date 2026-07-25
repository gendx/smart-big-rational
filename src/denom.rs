// Copyright 2026 The SmartBigRational Authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use num_bigint::{BigInt, BigUint};
use num_traits::{One, Pow};
use std::ops::{Div, DivAssign, Mul, MulAssign};

/// Interface representing the positive denominator of a rational number,
/// suitable for use in a [`SmartBigRational`](crate::SmartBigRational).
pub trait Denom:
    Clone
    + From<u8>
    + From<u16>
    + From<u32>
    + From<u64>
    + From<u128>
    + From<usize>
    + From<BigUint>
    + for<'a> From<&'a BigUint>
    + Into<BigUint>
    + One
    + Pow<u32, Output = Self>
    + Mul<Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + MulAssign
    + for<'a> MulAssign<&'a Self>
    + Div<Output = Self>
    + for<'a> Div<&'a Self, Output = Self>
    + DivAssign
    + for<'a> DivAssign<&'a Self>
    + Mul<BigInt, Output = BigInt>
    + for<'a> Mul<&'a BigInt, Output = BigInt>
{
    /// Constant value of 1.
    const ONE: Self;

    /// Converts this denominator into a big integer.
    fn into_biguint(self) -> BigUint;

    /// Converts this denominator into a big integer.
    fn to_biguint(&self) -> BigUint;

    /// Returns the least common multiple of two denominators, adjusting the
    /// numerators accordingly.
    fn normalize(lnum: &mut BigInt, rnum: &mut BigInt, ldenom: &Self, rdenom: &Self) -> Self;

    /// Reduces this denominator together with the given numerator so that their
    /// GCD is one.
    fn gcd_reduce(&mut self, num: &mut BigInt);
}

/// Additional trait that references to a [`Denom`] must implement.
pub trait DenomRef<D: Denom>:
    Into<BigUint>
    + Pow<u32, Output = D>
    + Mul<Self, Output = D>
    + Mul<D, Output = D>
    + Div<Self, Output = D>
    + Div<D, Output = D>
    + Mul<BigInt, Output = BigInt>
    + for<'a> Mul<&'a BigInt, Output = BigInt>
{
}
