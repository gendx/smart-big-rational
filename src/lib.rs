// Copyright 2023-2026 The SmartBigRational Authors
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

//! This crate implements a big rational type that is optimized for addition,
//! subtraction and multiplication.
//!
//! Unlike the vanilla [`BigRational`] type, the [`SmartBigRational`] doesn't
//! perform a full GCD reduction upon addition and subtraction, and doesn't
//! perform any reduction upon multiplication.

#![forbid(missing_docs, unsafe_code)]

mod denom;

pub use denom::Denom;
use num_bigint::{BigInt, BigUint, Sign};
use num_rational::BigRational;
use num_traits::{One, Pow, Zero};
use std::cmp::Ordering;
use std::fmt::Display;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// A big rational type optimized for addition, subtraction and multiplication.
///
/// This is achieved by representing the denominator with the [`Denom`] type,
/// and performing only partial GCD reductions during arithmetic operations.
#[derive(Clone, Debug)]
pub struct SmartBigRational {
    num: BigInt,
    denom: Denom,
}

impl SmartBigRational {
    /// Constant value of 0.
    pub const ZERO: Self = Self {
        num: BigInt::ZERO,
        denom: Denom::ONE,
    };

    /// Constant value of 1.
    pub const ONE: Self = Self {
        num: BigInt::ONE,
        denom: Denom::ONE,
    };

    /// Creates a new rational number by dividing the given numerator by the
    /// denominator.
    ///
    /// ```
    /// # use num_bigint::BigInt;
    /// # use num_rational::BigRational;
    /// # use smart_big_rational::SmartBigRational;
    /// let x = SmartBigRational::ratio(2, 3u32);
    ///
    /// assert_eq!(
    ///     BigRational::from(x),
    ///     BigRational::new(BigInt::from(2), BigInt::from(3))
    /// );
    /// ```
    pub fn ratio(num: impl Into<BigInt>, denom: impl Into<Denom>) -> Self {
        Self {
            num: num.into(),
            denom: denom.into(),
        }
    }

    /// Returns the current numerator and denominator as is, without reduction.
    pub fn into_raw(self) -> (BigInt, Denom) {
        (self.num, self.denom)
    }

    /// Returns the current numerator as is, without reduction.
    pub fn numer(&self) -> &BigInt {
        &self.num
    }

    /// Returns the current denominator as is, without reduction.
    pub fn denom(&self) -> &Denom {
        &self.denom
    }

    /// Converts this value to a [`BigRational`].
    pub fn into_big_rational(self) -> BigRational {
        self.into()
    }

    /// Converts this value to a [`BigRational`].
    pub fn to_big_rational(&self) -> BigRational {
        self.into()
    }

    /// Reduces the current value.
    ///
    /// After reduction, the GCD of the numerator and denominator is one.
    ///
    /// This is a slow operation, but may be beneficial in some cases (for
    /// example if this value is then used many times) as the representation
    /// becomes smaller if the numerator and denominator had many common
    /// factors. This may however be detrimental if you then add/subtract values
    /// that contain the same common factors that were reduced. Therefore,
    /// there is no rule of thumb: benchmark your concrete code to see if
    /// this brings any performance improvement.
    pub fn reduce(&mut self) {
        self.denom.gcd_reduce(&mut self.num);
    }
}

impl From<BigRational> for SmartBigRational {
    fn from(value: BigRational) -> SmartBigRational {
        let (num, denom) = value.into_raw();
        let (sign, denom) = denom.into_parts();
        assert_eq!(sign, Sign::Plus);
        SmartBigRational {
            num,
            denom: denom.into(),
        }
    }
}

impl From<&BigRational> for SmartBigRational {
    fn from(value: &BigRational) -> SmartBigRational {
        let denom = value.denom();
        assert_eq!(denom.sign(), Sign::Plus);
        SmartBigRational {
            num: value.numer().clone(),
            denom: denom.magnitude().into(),
        }
    }
}

impl From<BigInt> for SmartBigRational {
    fn from(value: BigInt) -> SmartBigRational {
        SmartBigRational {
            num: value,
            denom: Denom::ONE,
        }
    }
}

impl From<SmartBigRational> for BigRational {
    fn from(value: SmartBigRational) -> BigRational {
        BigRational::new(value.num, value.denom.to_biguint().into())
    }
}

impl From<&SmartBigRational> for BigRational {
    fn from(value: &SmartBigRational) -> BigRational {
        BigRational::new(value.num.clone(), value.denom.to_biguint().into())
    }
}

impl PartialEq for SmartBigRational {
    fn eq(&self, rhs: &Self) -> bool {
        self.num.sign() == rhs.num.sign()
            && self.num.magnitude() * rhs.denom.to_biguint()
                == rhs.num.magnitude() * self.denom.to_biguint()
    }
}

impl Eq for SmartBigRational {}

impl PartialOrd for SmartBigRational {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}

impl Ord for SmartBigRational {
    fn cmp(&self, rhs: &Self) -> Ordering {
        match (self.num.sign(), rhs.num.sign()) {
            (Sign::Plus, Sign::Plus) => (self.num.magnitude() * rhs.denom.to_biguint())
                .cmp(&(rhs.num.magnitude() * self.denom.to_biguint())),
            (Sign::Plus, Sign::NoSign) => Ordering::Greater,
            (Sign::Plus, Sign::Minus) => Ordering::Greater,
            (Sign::NoSign, Sign::Plus) => Ordering::Less,
            (Sign::NoSign, Sign::NoSign) => Ordering::Equal,
            (Sign::NoSign, Sign::Minus) => Ordering::Greater,
            (Sign::Minus, Sign::Plus) => Ordering::Less,
            (Sign::Minus, Sign::NoSign) => Ordering::Less,
            (Sign::Minus, Sign::Minus) => (rhs.num.magnitude() * self.denom.to_biguint())
                .cmp(&(self.num.magnitude() * rhs.denom.to_biguint())),
        }
    }
}

impl Display for SmartBigRational {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        Display::fmt(&self.to_big_rational(), f)
    }
}

impl Zero for SmartBigRational {
    fn zero() -> Self {
        Self::ZERO
    }

    fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
}

impl One for SmartBigRational {
    fn one() -> Self {
        Self::ONE
    }
}

impl Neg for SmartBigRational {
    type Output = Self;

    fn neg(self) -> Self {
        SmartBigRational {
            num: -self.num,
            denom: self.denom,
        }
    }
}

impl Neg for &SmartBigRational {
    type Output = SmartBigRational;

    fn neg(self) -> SmartBigRational {
        SmartBigRational {
            num: -&self.num,
            denom: self.denom.clone(),
        }
    }
}

impl Pow<u32> for SmartBigRational {
    type Output = Self;

    fn pow(self, rhs: u32) -> Self {
        SmartBigRational {
            num: self.num.pow(rhs),
            denom: self.denom.pow(rhs),
        }
    }
}

impl Pow<u32> for &SmartBigRational {
    type Output = SmartBigRational;

    fn pow(self, rhs: u32) -> SmartBigRational {
        SmartBigRational {
            num: Pow::pow(&self.num, rhs),
            denom: Pow::pow(&self.denom, rhs),
        }
    }
}

impl Add for SmartBigRational {
    type Output = Self;

    fn add(mut self, mut rhs: Self) -> Self {
        let denom = Denom::normalize(&mut self.num, &mut rhs.num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: self.num + rhs.num,
            denom,
        }
    }
}

impl Add<&SmartBigRational> for SmartBigRational {
    type Output = Self;

    fn add(mut self, rhs: &SmartBigRational) -> SmartBigRational {
        let mut rhs_num = rhs.num.clone();
        let denom = Denom::normalize(&mut self.num, &mut rhs_num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: self.num + rhs_num,
            denom,
        }
    }
}

impl Add for &SmartBigRational {
    type Output = SmartBigRational;

    fn add(self, rhs: Self) -> SmartBigRational {
        let mut num = self.num.clone();
        let mut rhs_num = rhs.num.clone();
        let denom = Denom::normalize(&mut num, &mut rhs_num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: num + rhs_num,
            denom,
        }
    }
}

impl Add<SmartBigRational> for &SmartBigRational {
    type Output = SmartBigRational;

    fn add(self, mut rhs: SmartBigRational) -> SmartBigRational {
        let mut num = self.num.clone();
        let denom = Denom::normalize(&mut num, &mut rhs.num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: num + rhs.num,
            denom,
        }
    }
}

impl AddAssign for SmartBigRational {
    fn add_assign(&mut self, mut rhs: Self) {
        self.denom = Denom::normalize(&mut self.num, &mut rhs.num, &self.denom, &rhs.denom);
        self.num += rhs.num;
    }
}

impl AddAssign<&SmartBigRational> for SmartBigRational {
    fn add_assign(&mut self, rhs: &SmartBigRational) {
        let mut rhs_num = rhs.num.clone();
        self.denom = Denom::normalize(&mut self.num, &mut rhs_num, &self.denom, &rhs.denom);
        self.num += rhs_num;
    }
}

impl Add<BigInt> for SmartBigRational {
    type Output = Self;

    fn add(self, mut rhs: BigInt) -> Self {
        rhs *= &self.denom;
        SmartBigRational {
            num: self.num + rhs,
            denom: self.denom,
        }
    }
}

impl Add<&BigInt> for SmartBigRational {
    type Output = Self;

    fn add(self, rhs: &BigInt) -> Self {
        SmartBigRational {
            num: self.num + rhs * &self.denom,
            denom: self.denom,
        }
    }
}

impl Add<BigInt> for &SmartBigRational {
    type Output = SmartBigRational;

    fn add(self, mut rhs: BigInt) -> SmartBigRational {
        rhs *= &self.denom;
        SmartBigRational {
            num: &self.num + rhs,
            denom: self.denom.clone(),
        }
    }
}

impl Add<&BigInt> for &SmartBigRational {
    type Output = SmartBigRational;

    fn add(self, rhs: &BigInt) -> SmartBigRational {
        SmartBigRational {
            num: &self.num + rhs * &self.denom,
            denom: self.denom.clone(),
        }
    }
}

impl AddAssign<BigInt> for SmartBigRational {
    fn add_assign(&mut self, mut rhs: BigInt) {
        rhs *= &self.denom;
        self.num += rhs;
    }
}

impl AddAssign<&BigInt> for SmartBigRational {
    fn add_assign(&mut self, rhs: &BigInt) {
        self.num += rhs * &self.denom;
    }
}

impl Sub for SmartBigRational {
    type Output = Self;

    fn sub(mut self, mut rhs: Self) -> Self {
        let denom = Denom::normalize(&mut self.num, &mut rhs.num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: self.num - rhs.num,
            denom,
        }
    }
}

impl Sub<&SmartBigRational> for SmartBigRational {
    type Output = Self;

    fn sub(mut self, rhs: &SmartBigRational) -> Self {
        let mut rhs_num = rhs.num.clone();
        let denom = Denom::normalize(&mut self.num, &mut rhs_num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: self.num - rhs_num,
            denom,
        }
    }
}

impl Sub for &SmartBigRational {
    type Output = SmartBigRational;

    fn sub(self, rhs: Self) -> SmartBigRational {
        let mut num = self.num.clone();
        let mut rhs_num = rhs.num.clone();
        let denom = Denom::normalize(&mut num, &mut rhs_num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: num - rhs_num,
            denom,
        }
    }
}

impl Sub<SmartBigRational> for &SmartBigRational {
    type Output = SmartBigRational;

    fn sub(self, mut rhs: SmartBigRational) -> SmartBigRational {
        let mut num = self.num.clone();
        let denom = Denom::normalize(&mut num, &mut rhs.num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: num - rhs.num,
            denom,
        }
    }
}

impl SubAssign for SmartBigRational {
    fn sub_assign(&mut self, mut rhs: Self) {
        self.denom = Denom::normalize(&mut self.num, &mut rhs.num, &self.denom, &rhs.denom);
        self.num -= rhs.num;
    }
}

impl SubAssign<&SmartBigRational> for SmartBigRational {
    fn sub_assign(&mut self, rhs: &SmartBigRational) {
        let mut rhs_num = rhs.num.clone();
        self.denom = Denom::normalize(&mut self.num, &mut rhs_num, &self.denom, &rhs.denom);
        self.num -= rhs_num;
    }
}

impl Sub<BigInt> for SmartBigRational {
    type Output = Self;

    fn sub(self, mut rhs: BigInt) -> Self {
        rhs *= &self.denom;
        SmartBigRational {
            num: self.num - rhs,
            denom: self.denom,
        }
    }
}

impl Sub<&BigInt> for SmartBigRational {
    type Output = Self;

    fn sub(self, rhs: &BigInt) -> Self {
        SmartBigRational {
            num: self.num - rhs * &self.denom,
            denom: self.denom,
        }
    }
}

impl Sub<BigInt> for &SmartBigRational {
    type Output = SmartBigRational;

    fn sub(self, mut rhs: BigInt) -> SmartBigRational {
        rhs *= &self.denom;
        SmartBigRational {
            num: &self.num - rhs,
            denom: self.denom.clone(),
        }
    }
}

impl Sub<&BigInt> for &SmartBigRational {
    type Output = SmartBigRational;

    fn sub(self, rhs: &BigInt) -> SmartBigRational {
        SmartBigRational {
            num: &self.num - rhs * &self.denom,
            denom: self.denom.clone(),
        }
    }
}

impl SubAssign<BigInt> for SmartBigRational {
    fn sub_assign(&mut self, mut rhs: BigInt) {
        rhs *= &self.denom;
        self.num -= rhs;
    }
}

impl SubAssign<&BigInt> for SmartBigRational {
    fn sub_assign(&mut self, rhs: &BigInt) {
        self.num -= rhs * &self.denom;
    }
}

impl Mul for SmartBigRational {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        SmartBigRational {
            num: self.num * rhs.num,
            denom: self.denom * rhs.denom,
        }
    }
}

impl Mul<&SmartBigRational> for SmartBigRational {
    type Output = Self;

    fn mul(self, rhs: &SmartBigRational) -> Self {
        SmartBigRational {
            num: self.num * &rhs.num,
            denom: self.denom * &rhs.denom,
        }
    }
}

impl Mul for &SmartBigRational {
    type Output = SmartBigRational;

    fn mul(self, rhs: Self) -> SmartBigRational {
        SmartBigRational {
            num: &self.num * &rhs.num,
            denom: &self.denom * &rhs.denom,
        }
    }
}

impl Mul<SmartBigRational> for &SmartBigRational {
    type Output = SmartBigRational;

    fn mul(self, rhs: SmartBigRational) -> SmartBigRational {
        SmartBigRational {
            num: &self.num * rhs.num,
            denom: &self.denom * rhs.denom,
        }
    }
}

impl MulAssign for SmartBigRational {
    fn mul_assign(&mut self, rhs: Self) {
        self.num *= rhs.num;
        self.denom *= rhs.denom;
    }
}

impl MulAssign<&SmartBigRational> for SmartBigRational {
    fn mul_assign(&mut self, rhs: &SmartBigRational) {
        self.num *= &rhs.num;
        self.denom *= &rhs.denom;
    }
}

impl Mul<BigInt> for SmartBigRational {
    type Output = Self;

    fn mul(self, rhs: BigInt) -> Self {
        SmartBigRational {
            num: self.num * rhs,
            denom: self.denom,
        }
    }
}

impl Mul<&BigInt> for SmartBigRational {
    type Output = Self;

    fn mul(self, rhs: &BigInt) -> Self {
        SmartBigRational {
            num: self.num * rhs,
            denom: self.denom,
        }
    }
}

impl Mul<BigInt> for &SmartBigRational {
    type Output = SmartBigRational;

    fn mul(self, rhs: BigInt) -> SmartBigRational {
        SmartBigRational {
            num: &self.num * rhs,
            denom: self.denom.clone(),
        }
    }
}

impl Mul<&BigInt> for &SmartBigRational {
    type Output = SmartBigRational;

    fn mul(self, rhs: &BigInt) -> SmartBigRational {
        SmartBigRational {
            num: &self.num * rhs,
            denom: self.denom.clone(),
        }
    }
}

impl MulAssign<BigInt> for SmartBigRational {
    fn mul_assign(&mut self, rhs: BigInt) {
        self.num *= rhs;
    }
}

impl MulAssign<&BigInt> for SmartBigRational {
    fn mul_assign(&mut self, rhs: &BigInt) {
        self.num *= rhs;
    }
}

impl Div for SmartBigRational {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        let (rhs_sign, rhs_num) = rhs.num.into_parts();
        let rhs_denom = BigInt::from_biguint(rhs_sign, BigUint::from(rhs.denom));
        SmartBigRational {
            num: self.num * rhs_denom,
            denom: self.denom * Denom::from(rhs_num),
        }
    }
}

impl Div<&SmartBigRational> for SmartBigRational {
    type Output = Self;

    fn div(self, rhs: &SmartBigRational) -> Self {
        let rhs_denom = BigInt::from_biguint(rhs.num.sign(), BigUint::from(&rhs.denom));
        SmartBigRational {
            num: self.num * rhs_denom,
            denom: self.denom * Denom::from(rhs.num.magnitude()),
        }
    }
}

impl Div for &SmartBigRational {
    type Output = SmartBigRational;

    fn div(self, rhs: Self) -> SmartBigRational {
        let rhs_denom = BigInt::from_biguint(rhs.num.sign(), BigUint::from(&rhs.denom));
        SmartBigRational {
            num: &self.num * rhs_denom,
            denom: &self.denom * Denom::from(rhs.num.magnitude()),
        }
    }
}

impl Div<SmartBigRational> for &SmartBigRational {
    type Output = SmartBigRational;

    fn div(self, rhs: SmartBigRational) -> SmartBigRational {
        let (rhs_sign, rhs_num) = rhs.num.into_parts();
        let rhs_denom = BigInt::from_biguint(rhs_sign, BigUint::from(rhs.denom));
        SmartBigRational {
            num: &self.num * rhs_denom,
            denom: &self.denom * Denom::from(rhs_num),
        }
    }
}

impl DivAssign for SmartBigRational {
    fn div_assign(&mut self, rhs: Self) {
        let (rhs_sign, rhs_num) = rhs.num.into_parts();
        let rhs_denom = BigInt::from_biguint(rhs_sign, BigUint::from(rhs.denom));
        self.num *= rhs_denom;
        self.denom *= Denom::from(rhs_num);
    }
}

impl DivAssign<&SmartBigRational> for SmartBigRational {
    fn div_assign(&mut self, rhs: &SmartBigRational) {
        let rhs_denom = BigInt::from_biguint(rhs.num.sign(), BigUint::from(&rhs.denom));
        self.num *= rhs_denom;
        self.denom *= Denom::from(rhs.num.magnitude());
    }
}

impl Div<BigUint> for SmartBigRational {
    type Output = Self;

    #[expect(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: BigUint) -> Self {
        SmartBigRational {
            num: self.num,
            denom: self.denom * Denom::from(rhs),
        }
    }
}

impl Div<&BigUint> for SmartBigRational {
    type Output = Self;

    #[expect(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: &BigUint) -> Self {
        SmartBigRational {
            num: self.num,
            denom: self.denom * Denom::from(rhs),
        }
    }
}

impl Div<BigUint> for &SmartBigRational {
    type Output = SmartBigRational;

    #[expect(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: BigUint) -> SmartBigRational {
        SmartBigRational {
            num: self.num.clone(),
            denom: &self.denom * Denom::from(rhs),
        }
    }
}

impl Div<&BigUint> for &SmartBigRational {
    type Output = SmartBigRational;

    #[expect(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: &BigUint) -> SmartBigRational {
        SmartBigRational {
            num: self.num.clone(),
            denom: &self.denom * Denom::from(rhs),
        }
    }
}

impl DivAssign<BigUint> for SmartBigRational {
    #[expect(clippy::suspicious_op_assign_impl)]
    fn div_assign(&mut self, rhs: BigUint) {
        self.denom *= Denom::from(rhs);
    }
}

impl DivAssign<&BigUint> for SmartBigRational {
    #[expect(clippy::suspicious_op_assign_impl)]
    fn div_assign(&mut self, rhs: &BigUint) {
        self.denom *= Denom::from(rhs);
    }
}

impl Sum for SmartBigRational {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a> Sum<&'a SmartBigRational> for SmartBigRational {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a SmartBigRational>,
    {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl Product for SmartBigRational {
    fn product<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a> Product<&'a SmartBigRational> for SmartBigRational {
    fn product<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a SmartBigRational>,
    {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::seq::IndexedRandom;

    fn get_positive_test_values() -> Vec<SmartBigRational> {
        let mut result = Vec::new();
        for i in 0..=30 {
            result.push(SmartBigRational::ratio(1 << i, 1u32));
        }
        for i in 0..=30 {
            result.push(SmartBigRational::ratio(1, 1u32 << i));
        }
        for i in 0..=30 {
            result.push(SmartBigRational::ratio(0x7FFF_FFFF - (1 << i), 1u32));
        }
        for i in 0..=30 {
            result.push(SmartBigRational::ratio(1, 0x7FFF_FFFF - (1u32 << i)));
        }
        result
    }

    fn loop_check1<T>(test_values: &[T], f: impl Fn(&T)) {
        for a in test_values {
            f(a);
        }
    }

    fn loop_check2<T>(test_values: &[T], f: impl Fn(&T, &T)) {
        for a in test_values {
            for b in test_values {
                f(a, b);
            }
        }
    }

    fn loop_check3<T>(test_values: &[T], num_samples: Option<usize>, f: impl Fn(&T, &T, &T)) {
        match num_samples {
            None => {
                // Exhaustive check.
                for a in test_values {
                    for b in test_values {
                        for c in test_values {
                            f(a, b, c);
                        }
                    }
                }
            }
            Some(n) => {
                // Randomly sample values rather than conducting an exhaustive O(n^3) search on
                // the test values.
                let mut rng = rand::rng();

                for _ in 0..n {
                    let a = test_values.choose(&mut rng).unwrap();
                    let b = test_values.choose(&mut rng).unwrap();
                    let c = test_values.choose(&mut rng).unwrap();
                    f(a, b, c);
                }
            }
        }
    }

    #[test]
    fn test_is_zero() {
        let test_values = get_positive_test_values();
        assert!(SmartBigRational::ZERO.is_zero());
        assert!(!SmartBigRational::ONE.is_zero());
        loop_check1(&test_values, |a| {
            assert!(!a.is_zero(), "{a} is zero");
        });
    }

    #[test]
    fn test_zero_is_add_neutral() {
        let test_values = get_positive_test_values();
        loop_check1(&test_values, |a| {
            assert_eq!(&(a + SmartBigRational::ZERO), a, "a + 0 != a for {a}");
            assert_eq!(&(SmartBigRational::ZERO + a), a, "0 + a != a for {a}");
            assert_eq!(&(a - SmartBigRational::ZERO), a, "a - 0 != a for {a}");
        })
    }

    #[test]
    fn test_add_is_commutative() {
        let test_values = get_positive_test_values();
        loop_check2(&test_values, |a, b| {
            assert_eq!(a + b, b + a, "a + b != b + a for {a}, {b}");
        })
    }

    #[test]
    fn test_add_is_associative() {
        let test_values = get_positive_test_values();
        loop_check3(&test_values, None, |a, b, c| {
            assert_eq!(
                (a + b) + c,
                a + (b + c),
                "(a + b) + c != a + (b + c) for {a}, {b}, {c}"
            );
        })
    }

    #[test]
    fn test_opposite() {
        let test_values = get_positive_test_values();
        loop_check1(&test_values, |a| {
            assert_eq!(&-(-a), a, "-(-a) != a for {a}");
            assert_eq!(
                &(SmartBigRational::ZERO - (SmartBigRational::ZERO - a)),
                a,
                "0 - (0 - a) != a for {a}"
            );
        });
    }

    #[test]
    fn test_sub_self() {
        let test_values = get_positive_test_values();
        loop_check1(&test_values, |a| {
            assert_eq!(a - a, SmartBigRational::ZERO, "a - a != 0 for {a}");
        });
    }

    #[test]
    fn test_add_sub() {
        let test_values = get_positive_test_values();
        loop_check2(&test_values, |a, b| {
            assert_eq!(&((a + b) - b), a, "(a + b) - b != a for {a}, {b}");
        });
    }

    #[test]
    fn test_sub_add() {
        let test_values = get_positive_test_values();
        loop_check2(&test_values, |a, b| {
            assert_eq!(&((a - b) + b), a, "(a - b) + b != a for {a}, {b}");
        });
    }

    #[test]
    fn test_one_is_mul_neutral() {
        let test_values = get_positive_test_values();
        loop_check1(&test_values, |a| {
            assert_eq!(&(a * SmartBigRational::ONE), a, "a * 1 != a for {a}");
            assert_eq!(&(SmartBigRational::ONE * a), a, "1 * a != a for {a}");
        })
    }

    #[test]
    fn test_mul_is_commutative() {
        let test_values = get_positive_test_values();
        loop_check2(&test_values, |a, b| {
            assert_eq!(a * b, b * a, "a * b != b * a for {a}, {b}");
        })
    }

    #[test]
    fn test_mul_is_associative() {
        let test_values = get_positive_test_values();
        loop_check3(&test_values, None, |a, b, c| {
            assert_eq!(
                (a * b) * c,
                a * (b * c),
                "(a * b) * c != a * (b * c) for {a}, {b}, {c}"
            );
        })
    }

    #[test]
    fn test_mul_is_distributive() {
        let test_values = get_positive_test_values();
        loop_check3(&test_values, None, |a, b, c| {
            assert_eq!(
                a * (b + c),
                (a * b) + (a * c),
                "a * (b + c) != (a * b) + (a * c) for {a}, {b}, {c}"
            );
        })
    }

    #[test]
    fn test_one_is_div_neutral() {
        let test_values = get_positive_test_values();
        loop_check1(&test_values, |a| {
            assert_eq!(&(a / SmartBigRational::ONE), a, "a / 1 != a for {a}");
        })
    }

    #[test]
    fn test_div_self() {
        let test_values = get_positive_test_values();
        loop_check1(&test_values, |a| {
            assert_eq!(a / a, SmartBigRational::ONE, "a / a != 1 for {a}");
        });
    }

    #[test]
    fn test_mul_div() {
        let test_values = get_positive_test_values();
        loop_check2(&test_values, |a, b| {
            assert_eq!(&((a * b) / b), a, "(a * b) / b != a for {a}, {b}");
        });
    }
}
