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
mod denom_array;
pub(crate) mod primes;

pub use denom::{Denom, DenomRef};
pub use denom_array::{DenomArray, DenomArray24};
use num_bigint::{BigInt, BigUint, Sign};
use num_rational::BigRational;
use num_traits::{One, Pow, Zero};
use std::cmp::Ordering;
use std::fmt::Display;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// A big rational type optimized for addition, subtraction and multiplication.
///
/// This is achieved by representing the denominator with a [`Denom`] type, and
/// performing only partial GCD reductions during arithmetic operations.
#[derive(Clone, Debug)]
pub struct SmartBigRational<D> {
    num: BigInt,
    denom: D,
}

impl<D: Denom> SmartBigRational<D> {
    /// Constant value of 0.
    pub const ZERO: Self = Self {
        num: BigInt::ZERO,
        denom: D::ONE,
    };

    /// Constant value of 1.
    pub const ONE: Self = Self {
        num: BigInt::ONE,
        denom: D::ONE,
    };

    /// Creates a new rational number by dividing the given numerator by the
    /// denominator.
    ///
    /// ```
    /// # use num_bigint::BigInt;
    /// # use num_rational::BigRational;
    /// # use smart_big_rational::{DenomArray24, SmartBigRational};
    /// let x = SmartBigRational::<DenomArray24>::ratio(2, 3u32);
    ///
    /// assert_eq!(
    ///     BigRational::from(x),
    ///     BigRational::new(BigInt::from(2), BigInt::from(3))
    /// );
    /// ```
    pub fn ratio(num: impl Into<BigInt>, denom: impl Into<D>) -> Self {
        Self {
            num: num.into(),
            denom: denom.into(),
        }
    }

    /// Returns the current numerator and denominator as is, without reduction.
    pub fn into_raw(self) -> (BigInt, D) {
        (self.num, self.denom)
    }

    /// Returns the current numerator as is, without reduction.
    pub fn numer(&self) -> &BigInt {
        &self.num
    }

    /// Returns the current denominator as is, without reduction.
    pub fn denom(&self) -> &D {
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

impl<D: Denom> From<BigRational> for SmartBigRational<D> {
    fn from(value: BigRational) -> Self {
        let (num, denom) = value.into_raw();
        let (sign, denom) = denom.into_parts();
        assert_eq!(sign, Sign::Plus);
        Self {
            num,
            denom: denom.into(),
        }
    }
}

impl<D: Denom> From<&BigRational> for SmartBigRational<D> {
    fn from(value: &BigRational) -> Self {
        let denom = value.denom();
        assert_eq!(denom.sign(), Sign::Plus);
        Self {
            num: value.numer().clone(),
            denom: denom.magnitude().into(),
        }
    }
}

impl<D: Denom> From<BigInt> for SmartBigRational<D> {
    fn from(value: BigInt) -> Self {
        Self {
            num: value,
            denom: D::ONE,
        }
    }
}

impl<D: Denom> From<SmartBigRational<D>> for BigRational {
    fn from(value: SmartBigRational<D>) -> BigRational {
        BigRational::new(value.num, value.denom.to_biguint().into())
    }
}

impl<D: Denom> From<&SmartBigRational<D>> for BigRational {
    fn from(value: &SmartBigRational<D>) -> BigRational {
        BigRational::new(value.num.clone(), value.denom.to_biguint().into())
    }
}

impl<D: Denom> PartialEq for SmartBigRational<D> {
    fn eq(&self, rhs: &Self) -> bool {
        self.num.sign() == rhs.num.sign()
            && self.num.magnitude() * rhs.denom.to_biguint()
                == rhs.num.magnitude() * self.denom.to_biguint()
    }
}

impl<D: Denom> Eq for SmartBigRational<D> {}

impl<D: Denom> PartialOrd for SmartBigRational<D> {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}

impl<D: Denom> Ord for SmartBigRational<D> {
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

impl<D: Denom> Display for SmartBigRational<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        Display::fmt(&self.to_big_rational(), f)
    }
}

impl<D: Denom> Zero for SmartBigRational<D> {
    fn zero() -> Self {
        Self::ZERO
    }

    fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
}

impl<D: Denom> One for SmartBigRational<D> {
    fn one() -> Self {
        Self::ONE
    }
}

impl<D: Denom> Neg for SmartBigRational<D> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            num: -self.num,
            denom: self.denom,
        }
    }
}

impl<D: Denom> Neg for &SmartBigRational<D> {
    type Output = SmartBigRational<D>;

    fn neg(self) -> SmartBigRational<D> {
        SmartBigRational {
            num: -&self.num,
            denom: self.denom.clone(),
        }
    }
}

impl<D: Denom> Pow<u32> for SmartBigRational<D> {
    type Output = Self;

    fn pow(self, rhs: u32) -> Self {
        Self {
            num: self.num.pow(rhs),
            denom: self.denom.pow(rhs),
        }
    }
}

impl<D: Denom> Pow<u32> for &SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = SmartBigRational<D>;

    fn pow(self, rhs: u32) -> SmartBigRational<D> {
        SmartBigRational {
            num: Pow::pow(&self.num, rhs),
            denom: Pow::pow(&self.denom, rhs),
        }
    }
}

impl<D: Denom> Add for SmartBigRational<D> {
    type Output = Self;

    fn add(mut self, mut rhs: Self) -> Self {
        let denom = D::normalize(&mut self.num, &mut rhs.num, &self.denom, &rhs.denom);
        Self {
            num: self.num + rhs.num,
            denom,
        }
    }
}

impl<D: Denom> Add<&Self> for SmartBigRational<D> {
    type Output = Self;

    fn add(mut self, rhs: &Self) -> Self {
        let mut rhs_num = rhs.num.clone();
        let denom = D::normalize(&mut self.num, &mut rhs_num, &self.denom, &rhs.denom);
        Self {
            num: self.num + rhs_num,
            denom,
        }
    }
}

impl<D: Denom> Add for &SmartBigRational<D> {
    type Output = SmartBigRational<D>;

    fn add(self, rhs: Self) -> SmartBigRational<D> {
        let mut num = self.num.clone();
        let mut rhs_num = rhs.num.clone();
        let denom = D::normalize(&mut num, &mut rhs_num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: num + rhs_num,
            denom,
        }
    }
}

impl<D: Denom> Add<SmartBigRational<D>> for &SmartBigRational<D> {
    type Output = SmartBigRational<D>;

    fn add(self, mut rhs: SmartBigRational<D>) -> SmartBigRational<D> {
        let mut num = self.num.clone();
        let denom = D::normalize(&mut num, &mut rhs.num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: num + rhs.num,
            denom,
        }
    }
}

impl<D: Denom> AddAssign for SmartBigRational<D> {
    fn add_assign(&mut self, mut rhs: Self) {
        self.denom = D::normalize(&mut self.num, &mut rhs.num, &self.denom, &rhs.denom);
        self.num += rhs.num;
    }
}

impl<D: Denom> AddAssign<&Self> for SmartBigRational<D> {
    fn add_assign(&mut self, rhs: &Self) {
        let mut rhs_num = rhs.num.clone();
        self.denom = D::normalize(&mut self.num, &mut rhs_num, &self.denom, &rhs.denom);
        self.num += rhs_num;
    }
}

impl<D: Denom> Add<BigInt> for SmartBigRational<D>
where
    BigInt: for<'a> MulAssign<&'a D>,
{
    type Output = Self;

    fn add(self, mut rhs: BigInt) -> Self {
        rhs *= &self.denom;
        Self {
            num: self.num + rhs,
            denom: self.denom,
        }
    }
}

impl<D: Denom> Add<&BigInt> for SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = Self;

    fn add(self, rhs: &BigInt) -> Self {
        Self {
            num: self.num + &self.denom * rhs,
            denom: self.denom,
        }
    }
}

impl<D: Denom> Add<BigInt> for &SmartBigRational<D>
where
    BigInt: for<'a> MulAssign<&'a D>,
{
    type Output = SmartBigRational<D>;

    fn add(self, mut rhs: BigInt) -> SmartBigRational<D> {
        rhs *= &self.denom;
        SmartBigRational {
            num: &self.num + rhs,
            denom: self.denom.clone(),
        }
    }
}

impl<D: Denom> Add<&BigInt> for &SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = SmartBigRational<D>;

    fn add(self, rhs: &BigInt) -> SmartBigRational<D> {
        SmartBigRational {
            num: &self.num + &self.denom * rhs,
            denom: self.denom.clone(),
        }
    }
}

impl<D: Denom> AddAssign<BigInt> for SmartBigRational<D>
where
    BigInt: for<'a> MulAssign<&'a D>,
{
    fn add_assign(&mut self, mut rhs: BigInt) {
        rhs *= &self.denom;
        self.num += rhs;
    }
}

impl<D: Denom> AddAssign<&BigInt> for SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    fn add_assign(&mut self, rhs: &BigInt) {
        self.num += &self.denom * rhs;
    }
}

impl<D: Denom> Sub for SmartBigRational<D> {
    type Output = Self;

    fn sub(mut self, mut rhs: Self) -> Self {
        let denom = D::normalize(&mut self.num, &mut rhs.num, &self.denom, &rhs.denom);
        Self {
            num: self.num - rhs.num,
            denom,
        }
    }
}

impl<D: Denom> Sub<&Self> for SmartBigRational<D> {
    type Output = Self;

    fn sub(mut self, rhs: &Self) -> Self {
        let mut rhs_num = rhs.num.clone();
        let denom = D::normalize(&mut self.num, &mut rhs_num, &self.denom, &rhs.denom);
        Self {
            num: self.num - rhs_num,
            denom,
        }
    }
}

impl<D: Denom> Sub for &SmartBigRational<D> {
    type Output = SmartBigRational<D>;

    fn sub(self, rhs: Self) -> SmartBigRational<D> {
        let mut num = self.num.clone();
        let mut rhs_num = rhs.num.clone();
        let denom = D::normalize(&mut num, &mut rhs_num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: num - rhs_num,
            denom,
        }
    }
}

impl<D: Denom> Sub<SmartBigRational<D>> for &SmartBigRational<D> {
    type Output = SmartBigRational<D>;

    fn sub(self, mut rhs: SmartBigRational<D>) -> SmartBigRational<D> {
        let mut num = self.num.clone();
        let denom = D::normalize(&mut num, &mut rhs.num, &self.denom, &rhs.denom);
        SmartBigRational {
            num: num - rhs.num,
            denom,
        }
    }
}

impl<D: Denom> SubAssign for SmartBigRational<D> {
    fn sub_assign(&mut self, mut rhs: Self) {
        self.denom = D::normalize(&mut self.num, &mut rhs.num, &self.denom, &rhs.denom);
        self.num -= rhs.num;
    }
}

impl<D: Denom> SubAssign<&Self> for SmartBigRational<D> {
    fn sub_assign(&mut self, rhs: &Self) {
        let mut rhs_num = rhs.num.clone();
        self.denom = D::normalize(&mut self.num, &mut rhs_num, &self.denom, &rhs.denom);
        self.num -= rhs_num;
    }
}

impl<D: Denom> Sub<BigInt> for SmartBigRational<D>
where
    BigInt: for<'a> MulAssign<&'a D>,
{
    type Output = Self;

    fn sub(self, mut rhs: BigInt) -> Self {
        rhs *= &self.denom;
        Self {
            num: self.num - rhs,
            denom: self.denom,
        }
    }
}

impl<D: Denom> Sub<&BigInt> for SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = Self;

    fn sub(self, rhs: &BigInt) -> Self {
        Self {
            num: self.num - &self.denom * rhs,
            denom: self.denom,
        }
    }
}

impl<D: Denom> Sub<BigInt> for &SmartBigRational<D>
where
    BigInt: for<'a> MulAssign<&'a D>,
{
    type Output = SmartBigRational<D>;

    fn sub(self, mut rhs: BigInt) -> SmartBigRational<D> {
        rhs *= &self.denom;
        SmartBigRational {
            num: &self.num - rhs,
            denom: self.denom.clone(),
        }
    }
}

impl<D: Denom> Sub<&BigInt> for &SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = SmartBigRational<D>;

    fn sub(self, rhs: &BigInt) -> SmartBigRational<D> {
        SmartBigRational {
            num: &self.num - &self.denom * rhs,
            denom: self.denom.clone(),
        }
    }
}

impl<D: Denom> SubAssign<BigInt> for SmartBigRational<D>
where
    BigInt: for<'a> MulAssign<&'a D>,
{
    fn sub_assign(&mut self, mut rhs: BigInt) {
        rhs *= &self.denom;
        self.num -= rhs;
    }
}

impl<D: Denom> SubAssign<&BigInt> for SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    fn sub_assign(&mut self, rhs: &BigInt) {
        self.num -= &self.denom * rhs;
    }
}

impl<D: Denom> Mul for SmartBigRational<D> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            num: self.num * rhs.num,
            denom: self.denom * rhs.denom,
        }
    }
}

impl<D: Denom> Mul<&Self> for SmartBigRational<D> {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self {
        Self {
            num: self.num * &rhs.num,
            denom: self.denom * &rhs.denom,
        }
    }
}

impl<D: Denom> Mul for &SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = SmartBigRational<D>;

    fn mul(self, rhs: Self) -> SmartBigRational<D> {
        SmartBigRational {
            num: &self.num * &rhs.num,
            denom: &self.denom * &rhs.denom,
        }
    }
}

impl<D: Denom> Mul<SmartBigRational<D>> for &SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = SmartBigRational<D>;

    fn mul(self, rhs: SmartBigRational<D>) -> SmartBigRational<D> {
        SmartBigRational {
            num: &self.num * rhs.num,
            denom: &self.denom * rhs.denom,
        }
    }
}

impl<D: Denom> MulAssign for SmartBigRational<D> {
    fn mul_assign(&mut self, rhs: Self) {
        self.num *= rhs.num;
        self.denom *= rhs.denom;
    }
}

impl<D: Denom> MulAssign<&Self> for SmartBigRational<D> {
    fn mul_assign(&mut self, rhs: &Self) {
        self.num *= &rhs.num;
        self.denom *= &rhs.denom;
    }
}

impl<D: Denom> Mul<BigInt> for SmartBigRational<D> {
    type Output = Self;

    fn mul(self, rhs: BigInt) -> Self {
        Self {
            num: self.num * rhs,
            denom: self.denom,
        }
    }
}

impl<D: Denom> Mul<&BigInt> for SmartBigRational<D> {
    type Output = Self;

    fn mul(self, rhs: &BigInt) -> Self {
        Self {
            num: self.num * rhs,
            denom: self.denom,
        }
    }
}

impl<D: Denom> Mul<BigInt> for &SmartBigRational<D> {
    type Output = SmartBigRational<D>;

    fn mul(self, rhs: BigInt) -> SmartBigRational<D> {
        SmartBigRational {
            num: &self.num * rhs,
            denom: self.denom.clone(),
        }
    }
}

impl<D: Denom> Mul<&BigInt> for &SmartBigRational<D> {
    type Output = SmartBigRational<D>;

    fn mul(self, rhs: &BigInt) -> SmartBigRational<D> {
        SmartBigRational {
            num: &self.num * rhs,
            denom: self.denom.clone(),
        }
    }
}

impl<D: Denom> MulAssign<BigInt> for SmartBigRational<D> {
    fn mul_assign(&mut self, rhs: BigInt) {
        self.num *= rhs;
    }
}

impl<D: Denom> MulAssign<&BigInt> for SmartBigRational<D> {
    fn mul_assign(&mut self, rhs: &BigInt) {
        self.num *= rhs;
    }
}

impl<D: Denom> Div for SmartBigRational<D> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        let (rhs_sign, rhs_num) = rhs.num.into_parts();
        let rhs_denom = BigInt::from_biguint(rhs_sign, rhs.denom.into());
        Self {
            num: self.num * rhs_denom,
            denom: self.denom * D::from(rhs_num),
        }
    }
}

impl<D: Denom> Div<&Self> for SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = Self;

    fn div(self, rhs: &Self) -> Self {
        let rhs_denom = BigInt::from_biguint(rhs.num.sign(), (&rhs.denom).into());
        Self {
            num: self.num * rhs_denom,
            denom: self.denom * D::from(rhs.num.magnitude()),
        }
    }
}

impl<D: Denom> Div for &SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = SmartBigRational<D>;

    fn div(self, rhs: Self) -> SmartBigRational<D> {
        let rhs_denom = BigInt::from_biguint(rhs.num.sign(), (&rhs.denom).into());
        SmartBigRational {
            num: &self.num * rhs_denom,
            denom: &self.denom * D::from(rhs.num.magnitude()),
        }
    }
}

impl<D: Denom> Div<SmartBigRational<D>> for &SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = SmartBigRational<D>;

    fn div(self, rhs: SmartBigRational<D>) -> SmartBigRational<D> {
        let (rhs_sign, rhs_num) = rhs.num.into_parts();
        let rhs_denom = BigInt::from_biguint(rhs_sign, rhs.denom.into());
        SmartBigRational {
            num: &self.num * rhs_denom,
            denom: &self.denom * D::from(rhs_num),
        }
    }
}

impl<D: Denom> DivAssign for SmartBigRational<D> {
    fn div_assign(&mut self, rhs: Self) {
        let (rhs_sign, rhs_num) = rhs.num.into_parts();
        let rhs_denom = BigInt::from_biguint(rhs_sign, rhs.denom.into());
        self.num *= rhs_denom;
        self.denom *= D::from(rhs_num);
    }
}

impl<D: Denom> DivAssign<&Self> for SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    fn div_assign(&mut self, rhs: &Self) {
        let rhs_denom = BigInt::from_biguint(rhs.num.sign(), (&rhs.denom).into());
        self.num *= rhs_denom;
        self.denom *= D::from(rhs.num.magnitude());
    }
}

impl<D: Denom> Div<BigUint> for SmartBigRational<D> {
    type Output = Self;

    #[expect(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: BigUint) -> Self {
        Self {
            num: self.num,
            denom: self.denom * D::from(rhs),
        }
    }
}

impl<D: Denom> Div<&BigUint> for SmartBigRational<D> {
    type Output = Self;

    #[expect(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: &BigUint) -> Self {
        Self {
            num: self.num,
            denom: self.denom * D::from(rhs),
        }
    }
}

impl<D: Denom> Div<BigUint> for &SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = SmartBigRational<D>;

    #[expect(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: BigUint) -> SmartBigRational<D> {
        SmartBigRational {
            num: self.num.clone(),
            denom: &self.denom * D::from(rhs),
        }
    }
}

impl<D: Denom> Div<&BigUint> for &SmartBigRational<D>
where
    for<'a> &'a D: DenomRef<D>,
{
    type Output = SmartBigRational<D>;

    #[expect(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: &BigUint) -> SmartBigRational<D> {
        SmartBigRational {
            num: self.num.clone(),
            denom: &self.denom * D::from(rhs),
        }
    }
}

impl<D: Denom> DivAssign<BigUint> for SmartBigRational<D> {
    #[expect(clippy::suspicious_op_assign_impl)]
    fn div_assign(&mut self, rhs: BigUint) {
        self.denom *= D::from(rhs);
    }
}

impl<D: Denom> DivAssign<&BigUint> for SmartBigRational<D> {
    #[expect(clippy::suspicious_op_assign_impl)]
    fn div_assign(&mut self, rhs: &BigUint) {
        self.denom *= D::from(rhs);
    }
}

impl<D: Denom> Sum for SmartBigRational<D> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<'a, D: Denom> Sum<&'a Self> for SmartBigRational<D> {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a Self>,
    {
        iter.fold(Self::zero(), |acc, x| acc + x)
    }
}

impl<D: Denom> Product for SmartBigRational<D> {
    fn product<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

impl<'a, D: Denom> Product<&'a Self> for SmartBigRational<D> {
    fn product<I>(iter: I) -> Self
    where
        I: Iterator<Item = &'a Self>,
    {
        iter.fold(Self::one(), |acc, x| acc * x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::seq::IndexedRandom;

    fn get_positive_test_values() -> Vec<SmartBigRational<DenomArray24>> {
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
        assert!(SmartBigRational::<DenomArray24>::ZERO.is_zero());
        assert!(!SmartBigRational::<DenomArray24>::ONE.is_zero());
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
