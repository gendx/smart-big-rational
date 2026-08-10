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

use crate::primes::ODD_PRIMES;
use crate::{Denom, DenomRef};
use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};
use std::cmp::Ordering;
use std::ops::{Div, DivAssign, Mul, MulAssign};

/// Denominator representation that decomposes an integer as a product of the
/// first `NUM_PRIMES` primes, multiplied by a regular big integer when that's
/// not sufficient.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DenomArray<const NUM_PRIMES: usize> {
    // Invariant: the representation is in canonical form, i.e. powers of small primes must
    // saturate the primes array before overflowing into the remainder.
    primes: [u8; NUM_PRIMES],
    remainder: Option<BigUint>,
}

/// Denominator representation that decomposes an integer as a product of the
/// first 24 primes (up to 89), multiplied by a regular big integer when that's
/// not sufficient.
pub type Denom24 = DenomArray<24>;

impl<const NUM_PRIMES: usize> DenomRef<DenomArray<NUM_PRIMES>> for &DenomArray<NUM_PRIMES> {}

impl<const NUM_PRIMES: usize> Denom for DenomArray<NUM_PRIMES> {
    const ONE: Self = Self {
        primes: [0; NUM_PRIMES],
        remainder: None,
    };

    fn into_biguint(self) -> BigUint {
        self.into()
    }

    fn to_biguint(&self) -> BigUint {
        self.into()
    }

    fn normalize(lnum: &mut BigInt, rnum: &mut BigInt, ldenom: &Self, rdenom: &Self) -> Self {
        let mut primes = [0; NUM_PRIMES];
        let mut ltmp = 1_usize;
        let mut rtmp = 1_usize;
        for i in 0..NUM_PRIMES {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            let lcount = ldenom.primes[i];
            let rcount = rdenom.primes[i];
            match lcount.cmp(&rcount) {
                Ordering::Equal => {
                    primes[i] = lcount;
                }
                Ordering::Less => {
                    Self::accum_pow(lnum, &mut ltmp, p, rcount - lcount);
                    primes[i] = rcount;
                }
                Ordering::Greater => {
                    Self::accum_pow(rnum, &mut rtmp, p, lcount - rcount);
                    primes[i] = lcount;
                }
            }
        }

        *lnum *= ltmp;
        *rnum *= rtmp;
        let remainder = match (&ldenom.remainder, &rdenom.remainder) {
            (None, None) => None,
            (None, Some(r)) => {
                *lnum *= BigInt::from_biguint(Sign::Plus, r.clone());
                Some(r.clone())
            }
            (Some(l), None) => {
                *rnum *= BigInt::from_biguint(Sign::Plus, l.clone());
                Some(l.clone())
            }
            (Some(l), Some(r)) => {
                if l == r {
                    Some(l.clone())
                } else {
                    *lnum *= BigInt::from_biguint(Sign::Plus, r.clone());
                    *rnum *= BigInt::from_biguint(Sign::Plus, l.clone());
                    Some(l * r)
                }
            }
        };
        Self { primes, remainder }
    }

    fn gcd_reduce(&mut self, num: &mut BigInt) {
        if let Some(remainder) = &mut self.remainder {
            let gcd = remainder.gcd(num.magnitude());
            if !gcd.is_one() {
                *remainder /= &gcd;
                if remainder.is_one() {
                    self.remainder = None;
                }
                let gcd: BigInt = gcd.into();
                *num /= gcd;
            }
        }

        'outer: for i in 0..NUM_PRIMES {
            let p = if i == 0 { 2 } else { ODD_PRIMES[i - 1] };
            if self.primes[i] != 0 {
                let p = BigInt::from(p);
                while self.primes[i] != 0 {
                    let (quo, rem) = num.div_rem(&p);
                    if !rem.is_zero() {
                        break;
                    }
                    *num = quo;
                    self.primes[i] -= 1;
                    if num.magnitude().is_one() {
                        break 'outer;
                    }
                }
            }
        }
    }
}

impl<const NUM_PRIMES: usize> DenomArray<NUM_PRIMES> {
    const _CHECK: () = assert!(NUM_PRIMES <= ODD_PRIMES.len());

    fn decompose_now(mut x: BigUint) -> Self {
        let mut primes = [0; NUM_PRIMES];

        let mut count2 = x.trailing_zeros().unwrap();
        if count2 != 0 {
            x >>= count2;
            if count2 <= u8::MAX as u64 {
                primes[0] = count2 as u8;
                count2 = 0;
            } else {
                primes[0] = u8::MAX;
                count2 -= u8::MAX as u64;
            }
        }

        'outer: for i in 1..NUM_PRIMES {
            let p = BigUint::from(ODD_PRIMES[i - 1]);
            while primes[i] != u8::MAX {
                let (quo, rem) = x.div_rem(&p);
                if !rem.is_zero() {
                    break;
                }
                x = quo;
                primes[i] += 1;
                if x.is_one() {
                    break 'outer;
                }
            }
        }

        let remainder = if x.is_one() && count2 == 0 {
            None
        } else {
            Some(x << count2)
        };
        Self { primes, remainder }
    }

    fn decompose_u8(mut x: u8) -> Self {
        let mut primes = [0; NUM_PRIMES];

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes[0] = count2 as u8;
        }

        'outer: for i in 1..NUM_PRIMES {
            let p = ODD_PRIMES[i - 1] as u8;
            while primes[i] != u8::MAX {
                if !x.is_multiple_of(p) {
                    break;
                }
                x /= p;
                primes[i] += 1;
                if x == 1 {
                    break 'outer;
                }
            }
        }

        let remainder = if x == 1 { None } else { Some(x.into()) };
        Self { primes, remainder }
    }

    fn decompose_u16(mut x: u16) -> Self {
        let mut primes = [0; NUM_PRIMES];

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes[0] = count2 as u8;
        }

        'outer: for i in 1..NUM_PRIMES {
            let p = ODD_PRIMES[i - 1];
            while primes[i] != u8::MAX {
                if !x.is_multiple_of(p) {
                    break;
                }
                x /= p;
                primes[i] += 1;
                if x == 1 {
                    break 'outer;
                }
            }
        }

        let remainder = if x == 1 { None } else { Some(x.into()) };
        Self { primes, remainder }
    }

    fn decompose_u32(mut x: u32) -> Self {
        let mut primes = [0; NUM_PRIMES];

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes[0] = count2 as u8;
        }

        'outer: for i in 1..NUM_PRIMES {
            let p = ODD_PRIMES[i - 1] as u32;
            while primes[i] != u8::MAX {
                if !x.is_multiple_of(p) {
                    break;
                }
                x /= p;
                primes[i] += 1;
                if x == 1 {
                    break 'outer;
                }
            }
        }

        let remainder = if x == 1 { None } else { Some(x.into()) };
        Self { primes, remainder }
    }

    fn decompose_u64(mut x: u64) -> Self {
        let mut primes = [0; NUM_PRIMES];

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes[0] = count2 as u8;
        }

        'outer: for i in 1..NUM_PRIMES {
            let p = ODD_PRIMES[i - 1] as u64;
            while primes[i] != u8::MAX {
                if !x.is_multiple_of(p) {
                    break;
                }
                x /= p;
                primes[i] += 1;
                if x == 1 {
                    break 'outer;
                }
            }
        }

        let remainder = if x == 1 { None } else { Some(x.into()) };
        Self { primes, remainder }
    }

    fn decompose_u128(mut x: u128) -> Self {
        let mut primes = [0; NUM_PRIMES];

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes[0] = count2 as u8;
        }

        'outer: for i in 1..NUM_PRIMES {
            let p = ODD_PRIMES[i - 1] as u128;
            while primes[i] != u8::MAX {
                if !x.is_multiple_of(p) {
                    break;
                }
                x /= p;
                primes[i] += 1;
                if x == 1 {
                    break 'outer;
                }
            }
        }

        let remainder = if x == 1 { None } else { Some(x.into()) };
        Self { primes, remainder }
    }

    fn decompose_usize(mut x: usize) -> Self {
        let mut primes = [0; NUM_PRIMES];

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes[0] = count2 as u8;
        }

        'outer: for i in 1..NUM_PRIMES {
            let p = ODD_PRIMES[i - 1] as usize;
            while primes[i] != u8::MAX {
                if !x.is_multiple_of(p) {
                    break;
                }
                x /= p;
                primes[i] += 1;
                if x == 1 {
                    break 'outer;
                }
            }
        }

        let remainder = if x == 1 { None } else { Some(x.into()) };
        Self { primes, remainder }
    }

    fn decompose_mask(
        remainder: &mut Option<BigUint>,
        primes: &mut [u8; NUM_PRIMES],
        mask: [bool; NUM_PRIMES],
    ) {
        let x: &mut BigUint = match remainder {
            None => return,
            Some(x) => x,
        };

        'outer: for i in 0..NUM_PRIMES {
            let p = if i == 0 { 2 } else { ODD_PRIMES[i - 1] };
            if !mask[i] || primes[i] == u8::MAX {
                continue;
            }
            let p = BigUint::from(p);
            while primes[i] != u8::MAX {
                let (quo, rem) = x.div_rem(&p);
                if !rem.is_zero() {
                    break;
                }
                *x = quo;
                primes[i] += 1;
                if x.is_one() {
                    break 'outer;
                }
            }
        }

        if x.is_one() {
            *remainder = None;
        }
    }

    /// Computes `prime.pow(exponent)` and multiplies it into the accumulated
    /// `(numerator, tmp)`.
    fn accum_pow(numerator: &mut BigInt, tmp: &mut usize, prime: usize, exponent: u8) {
        for _ in 0..exponent {
            match tmp.checked_mul(prime) {
                Some(prod) => *tmp = prod,
                None => {
                    *numerator *= *tmp;
                    *tmp = prime;
                }
            }
        }
    }

    fn pow_primes(
        this: &[u8; NUM_PRIMES],
        exponent: u32,
        remainder: &mut Option<BigUint>,
    ) -> [u8; NUM_PRIMES] {
        let mut primes = [0; NUM_PRIMES];
        for i in 0..NUM_PRIMES {
            let p = if i == 0 { 2 } else { ODD_PRIMES[i - 1] };
            let product = (this[i] as u32).strict_mul(exponent);
            if product <= u8::MAX as u32 {
                primes[i] = product as u8;
            } else {
                primes[i] = u8::MAX;
                let factor = BigUint::from(p).pow(product - u8::MAX as u32);
                match remainder {
                    None => *remainder = Some(factor),
                    Some(r) => *r *= factor,
                };
            }
        }
        primes
    }

    fn mul_primes(
        lhs: &[u8; NUM_PRIMES],
        rhs: &[u8; NUM_PRIMES],
        remainder: &mut Option<BigUint>,
    ) -> [u8; NUM_PRIMES] {
        let mut primes = [0; NUM_PRIMES];
        for i in 0..NUM_PRIMES {
            let p = if i == 0 { 2 } else { ODD_PRIMES[i - 1] };
            let sum = lhs[i] as u32 + rhs[i] as u32;
            if sum <= u8::MAX as u32 {
                primes[i] = sum as u8;
            } else {
                primes[i] = u8::MAX;
                let factor = BigUint::from(p).pow(sum - u8::MAX as u32);
                match remainder {
                    None => *remainder = Some(factor),
                    Some(r) => *r *= factor,
                };
            }
        }
        primes
    }

    fn div_primes(
        num: &[u8; NUM_PRIMES],
        denom: &[u8; NUM_PRIMES],
        mut remainder: Option<BigUint>,
    ) -> (Self, [bool; NUM_PRIMES]) {
        let mut primes = [0; NUM_PRIMES];
        let mut mask = [false; NUM_PRIMES];
        for i in 0..NUM_PRIMES {
            let p = if i == 0 { 2 } else { ODD_PRIMES[i - 1] };
            mask[i] = num[i] == u8::MAX;
            if num[i] >= denom[i] {
                primes[i] = num[i] - denom[i];
            } else {
                primes[i] = 0;
                let factor = BigUint::from(p).pow(denom[i] as u32 - num[i] as u32);
                remainder = match remainder {
                    None => Some(factor),
                    Some(r) => Some(r * factor),
                };
            }
        }
        (Self { primes, remainder }, mask)
    }
}

impl<const NUM_PRIMES: usize> From<u8> for DenomArray<NUM_PRIMES> {
    fn from(value: u8) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_u8(value)
    }
}

impl<const NUM_PRIMES: usize> From<u16> for DenomArray<NUM_PRIMES> {
    fn from(value: u16) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_u16(value)
    }
}

impl<const NUM_PRIMES: usize> From<u32> for DenomArray<NUM_PRIMES> {
    fn from(value: u32) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_u32(value)
    }
}

impl<const NUM_PRIMES: usize> From<u64> for DenomArray<NUM_PRIMES> {
    fn from(value: u64) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_u64(value)
    }
}

impl<const NUM_PRIMES: usize> From<u128> for DenomArray<NUM_PRIMES> {
    fn from(value: u128) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_u128(value)
    }
}

impl<const NUM_PRIMES: usize> From<usize> for DenomArray<NUM_PRIMES> {
    fn from(value: usize) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_usize(value)
    }
}

impl<const NUM_PRIMES: usize> From<BigUint> for DenomArray<NUM_PRIMES> {
    fn from(value: BigUint) -> Self {
        if value.is_zero() {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_now(value)
    }
}

impl<const NUM_PRIMES: usize> From<&BigUint> for DenomArray<NUM_PRIMES> {
    fn from(value: &BigUint) -> Self {
        if value.is_zero() {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_now(value.clone())
    }
}

impl<const NUM_PRIMES: usize> From<DenomArray<NUM_PRIMES>> for BigUint {
    fn from(value: DenomArray<NUM_PRIMES>) -> BigUint {
        let mut result = match value.remainder {
            Some(x) => x,
            None => BigUint::ONE,
        };
        let mut tmp = 1_usize;
        for (i, &count) in value.primes.iter().enumerate() {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            for _ in 0..count {
                match tmp.checked_mul(p) {
                    Some(prod) => tmp = prod,
                    None => {
                        result *= tmp;
                        tmp = p;
                    }
                }
            }
        }
        result * tmp
    }
}

impl<const NUM_PRIMES: usize> From<&DenomArray<NUM_PRIMES>> for BigUint {
    fn from(value: &DenomArray<NUM_PRIMES>) -> BigUint {
        let mut result = match &value.remainder {
            Some(x) => x.clone(),
            None => BigUint::ONE,
        };
        let mut tmp = 1_usize;
        for (i, &count) in value.primes.iter().enumerate() {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            for _ in 0..count {
                match tmp.checked_mul(p) {
                    Some(prod) => tmp = prod,
                    None => {
                        result *= tmp;
                        tmp = p;
                    }
                }
            }
        }
        result * tmp
    }
}

impl<const NUM_PRIMES: usize> One for DenomArray<NUM_PRIMES> {
    fn one() -> Self {
        Self::ONE
    }
}

impl<const NUM_PRIMES: usize> num_traits::Pow<u32> for DenomArray<NUM_PRIMES> {
    type Output = Self;

    fn pow(mut self, rhs: u32) -> Self {
        if rhs == 0 {
            return Self::ONE;
        }
        let primes = Self::pow_primes(&self.primes, rhs, &mut self.remainder);
        Self {
            primes,
            remainder: self.remainder,
        }
    }
}

impl<const NUM_PRIMES: usize> num_traits::Pow<u32> for &DenomArray<NUM_PRIMES> {
    type Output = DenomArray<NUM_PRIMES>;

    fn pow(self, rhs: u32) -> DenomArray<NUM_PRIMES> {
        if rhs == 0 {
            return DenomArray::ONE;
        }
        let mut remainder = self.remainder.clone();
        let primes = DenomArray::pow_primes(&self.primes, rhs, &mut remainder);
        DenomArray { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Mul for DenomArray<NUM_PRIMES> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut remainder = match (self.remainder, rhs.remainder) {
            (None, None) => None,
            (None, Some(r)) => Some(r),
            (Some(l), None) => Some(l),
            (Some(l), Some(r)) => Some(l * r),
        };
        let primes = Self::mul_primes(&self.primes, &rhs.primes, &mut remainder);
        Self { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Mul<&Self> for DenomArray<NUM_PRIMES> {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self {
        let mut remainder = match (self.remainder, &rhs.remainder) {
            (None, None) => None,
            (None, Some(r)) => Some(r.clone()),
            (Some(l), None) => Some(l),
            (Some(l), Some(r)) => Some(l * r),
        };
        let primes = Self::mul_primes(&self.primes, &rhs.primes, &mut remainder);
        Self { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Mul for &DenomArray<NUM_PRIMES> {
    type Output = DenomArray<NUM_PRIMES>;

    fn mul(self, rhs: Self) -> DenomArray<NUM_PRIMES> {
        let mut remainder = match (&self.remainder, &rhs.remainder) {
            (None, None) => None,
            (None, Some(r)) => Some(r.clone()),
            (Some(l), None) => Some(l.clone()),
            (Some(l), Some(r)) => Some(l * r),
        };
        let primes = DenomArray::mul_primes(&self.primes, &rhs.primes, &mut remainder);
        DenomArray { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Mul<DenomArray<NUM_PRIMES>> for &DenomArray<NUM_PRIMES> {
    type Output = DenomArray<NUM_PRIMES>;

    fn mul(self, rhs: DenomArray<NUM_PRIMES>) -> DenomArray<NUM_PRIMES> {
        let mut remainder = match (&self.remainder, rhs.remainder) {
            (None, None) => None,
            (None, Some(r)) => Some(r),
            (Some(l), None) => Some(l.clone()),
            (Some(l), Some(r)) => Some(l * r),
        };
        let primes = DenomArray::mul_primes(&self.primes, &rhs.primes, &mut remainder);
        DenomArray { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> MulAssign for DenomArray<NUM_PRIMES> {
    fn mul_assign(&mut self, rhs: Self) {
        match (&mut self.remainder, rhs.remainder) {
            (_, None) => (),
            (None, Some(r)) => self.remainder = Some(r),
            (Some(l), Some(r)) => *l *= r,
        };
        self.primes = Self::mul_primes(&self.primes, &rhs.primes, &mut self.remainder);
    }
}

impl<const NUM_PRIMES: usize> MulAssign<&Self> for DenomArray<NUM_PRIMES> {
    fn mul_assign(&mut self, rhs: &Self) {
        match (&mut self.remainder, &rhs.remainder) {
            (_, None) => (),
            (None, Some(r)) => self.remainder = Some(r.clone()),
            (Some(l), Some(r)) => *l *= r,
        };
        self.primes = Self::mul_primes(&self.primes, &rhs.primes, &mut self.remainder);
    }
}

impl<const NUM_PRIMES: usize> Div for DenomArray<NUM_PRIMES> {
    type Output = Self;

    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div(self, rhs: Self) -> Self {
        let (
            Self {
                mut primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = Self::div_primes(&self.primes, &rhs.primes, rhs.remainder);

        let mut remainder = match (self.remainder, rhs_remainder) {
            (None, None) => None,
            (None, Some(r)) => {
                if !r.is_one() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
                None
            }
            (Some(l), None) => Some(l),
            (Some(l), Some(r)) => {
                let (quo, rem) = l.div_rem(&r);
                if !rem.is_zero() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
                if quo.is_one() { None } else { Some(quo) }
            }
        };

        Self::decompose_mask(&mut remainder, &mut primes, mask);
        Self { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Div<&Self> for DenomArray<NUM_PRIMES> {
    type Output = Self;

    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div(self, rhs: &Self) -> Self {
        let (
            Self {
                mut primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = Self::div_primes(&self.primes, &rhs.primes, rhs.remainder.clone());

        let mut remainder = match (self.remainder, rhs_remainder) {
            (None, None) => None,
            (None, Some(r)) => {
                if !r.is_one() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
                None
            }
            (Some(l), None) => Some(l),
            (Some(l), Some(r)) => {
                let (quo, rem) = l.div_rem(&r);
                if !rem.is_zero() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
                if quo.is_one() { None } else { Some(quo) }
            }
        };

        Self::decompose_mask(&mut remainder, &mut primes, mask);
        Self { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Div for &DenomArray<NUM_PRIMES> {
    type Output = DenomArray<NUM_PRIMES>;

    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div(self, rhs: Self) -> DenomArray<NUM_PRIMES> {
        let (
            DenomArray {
                mut primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = DenomArray::div_primes(&self.primes, &rhs.primes, rhs.remainder.clone());

        let mut remainder = match (&self.remainder, rhs_remainder) {
            (None, None) => None,
            (None, Some(r)) => {
                if !r.is_one() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
                None
            }
            (Some(l), None) => Some(l.clone()),
            (Some(l), Some(r)) => {
                let (quo, rem) = l.div_rem(&r);
                if !rem.is_zero() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
                if quo.is_one() { None } else { Some(quo) }
            }
        };

        DenomArray::decompose_mask(&mut remainder, &mut primes, mask);
        DenomArray { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Div<DenomArray<NUM_PRIMES>> for &DenomArray<NUM_PRIMES> {
    type Output = DenomArray<NUM_PRIMES>;

    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div(self, rhs: DenomArray<NUM_PRIMES>) -> DenomArray<NUM_PRIMES> {
        let (
            DenomArray {
                mut primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = DenomArray::div_primes(&self.primes, &rhs.primes, rhs.remainder);

        let mut remainder = match (&self.remainder, rhs_remainder) {
            (None, None) => None,
            (None, Some(r)) => {
                if !r.is_one() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
                None
            }
            (Some(l), None) => Some(l.clone()),
            (Some(l), Some(r)) => {
                let (quo, rem) = l.div_rem(&r);
                if !rem.is_zero() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
                if quo.is_one() { None } else { Some(quo) }
            }
        };

        DenomArray::decompose_mask(&mut remainder, &mut primes, mask);
        DenomArray { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> DivAssign for DenomArray<NUM_PRIMES> {
    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div_assign(&mut self, rhs: Self) {
        let (
            Self {
                primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = Self::div_primes(&self.primes, &rhs.primes, rhs.remainder);
        self.primes = primes;

        match (&mut self.remainder, rhs_remainder) {
            (_, None) => (),
            (None, Some(r)) => {
                if !r.is_one() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
            }
            (Some(l), Some(r)) => {
                let (quo, rem) = l.div_rem(&r);
                if !rem.is_zero() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
                if quo.is_one() {
                    self.remainder = None
                } else {
                    *l = quo;
                }
            }
        };

        Self::decompose_mask(&mut self.remainder, &mut self.primes, mask);
    }
}

impl<const NUM_PRIMES: usize> DivAssign<&Self> for DenomArray<NUM_PRIMES> {
    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div_assign(&mut self, rhs: &Self) {
        let (
            Self {
                primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = Self::div_primes(&self.primes, &rhs.primes, rhs.remainder.clone());
        self.primes = primes;

        match (&mut self.remainder, rhs_remainder) {
            (_, None) => (),
            (None, Some(r)) => {
                if !r.is_one() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
            }
            (Some(l), Some(r)) => {
                let (quo, rem) = l.div_rem(&r);
                if !rem.is_zero() {
                    panic!("Attempted to divide a denominator by a non-divisor");
                }
                if quo.is_one() {
                    self.remainder = None
                } else {
                    *l = quo;
                }
            }
        };

        Self::decompose_mask(&mut self.remainder, &mut self.primes, mask);
    }
}

impl<const NUM_PRIMES: usize> Mul<DenomArray<NUM_PRIMES>> for BigInt {
    type Output = Self;

    fn mul(mut self, rhs: DenomArray<NUM_PRIMES>) -> Self {
        let mut tmp = 1_usize;
        for i in 0..NUM_PRIMES {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            let count = rhs.primes[i];
            if count != 0 {
                DenomArray::<NUM_PRIMES>::accum_pow(&mut self, &mut tmp, p, count);
            }
        }

        self *= tmp;
        if let Some(remainder) = rhs.remainder {
            self *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        self
    }
}

impl<const NUM_PRIMES: usize> Mul<&DenomArray<NUM_PRIMES>> for BigInt {
    type Output = Self;

    fn mul(mut self, rhs: &DenomArray<NUM_PRIMES>) -> Self {
        let mut tmp = 1_usize;
        for i in 0..NUM_PRIMES {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            let count = rhs.primes[i];
            if count != 0 {
                DenomArray::<NUM_PRIMES>::accum_pow(&mut self, &mut tmp, p, count);
            }
        }

        self *= tmp;
        if let Some(remainder) = &rhs.remainder {
            self *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        self
    }
}

impl<const NUM_PRIMES: usize> Mul<DenomArray<NUM_PRIMES>> for &BigInt {
    type Output = BigInt;

    fn mul(self, rhs: DenomArray<NUM_PRIMES>) -> BigInt {
        let mut this = self.clone();
        let mut tmp = 1_usize;
        for i in 0..NUM_PRIMES {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            let count = rhs.primes[i];
            if count != 0 {
                DenomArray::<NUM_PRIMES>::accum_pow(&mut this, &mut tmp, p, count);
            }
        }

        this *= tmp;
        if let Some(remainder) = rhs.remainder {
            this *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        this
    }
}

impl<const NUM_PRIMES: usize> Mul<&DenomArray<NUM_PRIMES>> for &BigInt {
    type Output = BigInt;

    fn mul(self, rhs: &DenomArray<NUM_PRIMES>) -> BigInt {
        let mut this = self.clone();
        let mut tmp = 1_usize;
        for i in 0..NUM_PRIMES {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            let count = rhs.primes[i];
            if count != 0 {
                DenomArray::<NUM_PRIMES>::accum_pow(&mut this, &mut tmp, p, count);
            }
        }

        this *= tmp;
        if let Some(remainder) = &rhs.remainder {
            this *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        this
    }
}

impl<const NUM_PRIMES: usize> Mul<BigInt> for DenomArray<NUM_PRIMES> {
    type Output = BigInt;

    fn mul(self, mut rhs: BigInt) -> BigInt {
        let mut tmp = 1_usize;
        for i in 0..NUM_PRIMES {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            let count = self.primes[i];
            if count != 0 {
                Self::accum_pow(&mut rhs, &mut tmp, p, count);
            }
        }

        rhs *= tmp;
        if let Some(remainder) = self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        rhs
    }
}

impl<const NUM_PRIMES: usize> Mul<&BigInt> for DenomArray<NUM_PRIMES> {
    type Output = BigInt;

    fn mul(self, rhs: &BigInt) -> BigInt {
        let mut rhs = rhs.clone();
        let mut tmp = 1_usize;
        for i in 0..NUM_PRIMES {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            let count = self.primes[i];
            if count != 0 {
                Self::accum_pow(&mut rhs, &mut tmp, p, count);
            }
        }

        rhs *= tmp;
        if let Some(remainder) = self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        rhs
    }
}

impl<const NUM_PRIMES: usize> Mul<BigInt> for &DenomArray<NUM_PRIMES> {
    type Output = BigInt;

    fn mul(self, mut rhs: BigInt) -> BigInt {
        let mut tmp = 1_usize;
        for i in 0..NUM_PRIMES {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            let count = self.primes[i];
            if count != 0 {
                DenomArray::<NUM_PRIMES>::accum_pow(&mut rhs, &mut tmp, p, count);
            }
        }

        rhs *= tmp;
        if let Some(remainder) = &self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        rhs
    }
}

impl<const NUM_PRIMES: usize> Mul<&BigInt> for &DenomArray<NUM_PRIMES> {
    type Output = BigInt;

    fn mul(self, rhs: &BigInt) -> BigInt {
        let mut rhs = rhs.clone();
        let mut tmp = 1_usize;
        for i in 0..NUM_PRIMES {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            let count = self.primes[i];
            if count != 0 {
                DenomArray::<NUM_PRIMES>::accum_pow(&mut rhs, &mut tmp, p, count);
            }
        }

        rhs *= tmp;
        if let Some(remainder) = &self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        rhs
    }
}

impl<const NUM_PRIMES: usize> MulAssign<DenomArray<NUM_PRIMES>> for BigInt {
    fn mul_assign(&mut self, rhs: DenomArray<NUM_PRIMES>) {
        let mut tmp = 1_usize;
        for i in 0..NUM_PRIMES {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            let count = rhs.primes[i];
            if count != 0 {
                DenomArray::<NUM_PRIMES>::accum_pow(self, &mut tmp, p, count);
            }
        }

        *self *= tmp;
        if let Some(remainder) = rhs.remainder {
            *self *= BigInt::from_biguint(Sign::Plus, remainder);
        }
    }
}

impl<const NUM_PRIMES: usize> MulAssign<&DenomArray<NUM_PRIMES>> for BigInt {
    fn mul_assign(&mut self, rhs: &DenomArray<NUM_PRIMES>) {
        let mut tmp = 1_usize;
        for i in 0..NUM_PRIMES {
            let p = if i == 0 {
                2
            } else {
                ODD_PRIMES[i - 1] as usize
            };
            let count = rhs.primes[i];
            if count != 0 {
                DenomArray::<NUM_PRIMES>::accum_pow(self, &mut tmp, p, count);
            }
        }

        *self *= tmp;
        if let Some(remainder) = &rhs.remainder {
            *self *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decompose_is_correct() {
        for i in 1_usize..=1000 {
            let bigi = BigUint::from(i);
            let x = Denom24::from(&bigi);
            let mut recomposed = x.remainder.unwrap_or_else(BigUint::one);
            for i in 0..24 {
                let prime = if i == 0 { 2 } else { ODD_PRIMES[i - 1] };
                for _ in 0..x.primes[i] {
                    recomposed *= prime;
                }
            }
            assert_eq!(recomposed, bigi);
        }
    }

    #[test]
    fn test_decompose_known_values() {
        assert_eq!(
            Denom24::from(BigUint::from(128_usize)),
            Denom24 {
                primes: [
                    7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                ],
                remainder: None,
            }
        );
        assert_eq!(
            Denom24::from(BigUint::from(89_usize)),
            Denom24 {
                primes: [
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
                ],
                remainder: None,
            }
        );
        assert_eq!(
            Denom24::from(BigUint::from(97_usize)),
            Denom24 {
                primes: [
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                ],
                remainder: Some(BigUint::from(97_usize)),
            }
        );
        assert_eq!(
            Denom24::from(BigUint::from(97000_usize)),
            Denom24 {
                primes: [
                    3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                ],
                remainder: Some(BigUint::from(97_usize)),
            }
        );
    }

    #[test]
    fn test_decompose_prime_powers() {
        for i in 0..24 {
            let p = if i == 0 { 2 } else { ODD_PRIMES[i - 1] };
            let p = BigUint::from(p);
            for power in 1..=255 {
                assert_eq!(
                    Denom24::from(p.pow(power as u32)),
                    Denom24 {
                        primes: std::array::from_fn(|j| if i == j { power } else { 0 }),
                        remainder: None,
                    }
                );
            }
            for power in 1..=64 {
                assert_eq!(
                    Denom24::from(p.pow(255 + power)),
                    Denom24 {
                        primes: std::array::from_fn(|j| if i == j { 255 } else { 0 }),
                        remainder: Some(p.pow(power)),
                    }
                );
            }
        }
    }

    #[test]
    fn test_mul_prime_powers() {
        let p = BigUint::from(2u32);
        for a in 1..=256 {
            let denom_a = Denom24::from(p.pow(a));
            for b in 1..=256 {
                let denom_b = Denom24::from(p.pow(b));
                let denom_ab = Denom24::from(p.pow(a + b));
                assert_eq!(&denom_a * denom_b, denom_ab);
            }
        }
    }

    #[test]
    fn test_div_prime_powers() {
        let p = BigUint::from(2u32);
        for a in 1..=256 {
            let denom_a = Denom24::from(p.pow(a));
            for b in 1..=256 {
                let denom_b = Denom24::from(p.pow(b));
                let denom_ab = Denom24::from(p.pow(a + b));
                assert_eq!(denom_ab / denom_b, denom_a);
            }
        }
    }

    #[test]
    fn test_to_biguint() {
        for i in 1_usize..=1000 {
            let bigi = BigUint::from(i);
            let x = Denom24::from(&bigi);
            assert_eq!(x.to_biguint(), bigi);
        }
    }

    #[test]
    fn test_product() {
        let values = (100..200)
            .map(|i: usize| Denom24::from(BigUint::from(i)))
            .collect::<Vec<_>>();
        for (i, x) in values.iter().enumerate().map(|(i, x)| (i + 100, x)) {
            for (j, y) in values.iter().enumerate().map(|(j, y)| (j + 100, y)) {
                let z = x * y;
                assert_eq!(z, Denom24::from(BigUint::from(i * j)));
                for k in 0..24 {
                    assert_eq!(z.primes[k], x.primes[k] + y.primes[k]);
                }
            }
        }
    }

    #[test]
    fn test_normalize() {
        let values = (100..200)
            .map(|i: usize| Denom24::from(BigUint::from(i)))
            .collect::<Vec<_>>();
        for x in &values {
            for y in &values {
                let mut xnum = BigInt::one();
                let mut ynum = BigInt::one();
                let lcm = Denom24::normalize(&mut xnum, &mut ynum, x, y);
                let lcm_bigint = lcm.to_biguint();
                let xnum = xnum.to_biguint().unwrap();
                let ynum = ynum.to_biguint().unwrap();

                assert_eq!(xnum * x.to_biguint(), lcm_bigint);
                assert_eq!(ynum * y.to_biguint(), lcm_bigint);
                for k in 0..24 {
                    assert_eq!(lcm.primes[k], std::cmp::max(x.primes[k], y.primes[k]));
                }
            }
        }
    }

    #[test]
    fn test_gcd_reduce() {
        let mut num = BigInt::from(-3 * 97);
        let mut denom = Denom24::from(BigUint::from(3u32 * 5 * 97));
        denom.gcd_reduce(&mut num);
        assert_eq!(num, BigInt::from(-1));
        assert_eq!(denom, Denom24::from(BigUint::from(5u32)));
    }
}
