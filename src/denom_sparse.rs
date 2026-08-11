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

use crate::primes::{
    ODD_PRIME_DIVIDERS_U16, ODD_PRIME_DIVIDERS_U32, ODD_PRIME_DIVIDERS_U64, ODD_PRIMES,
    known_odd_prime_factor_indices,
};
use crate::util::OddDivider;
use crate::{Denom, DenomRef};
use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};
use std::cmp::Ordering;
use std::iter::Peekable;
use std::ops::{Div, DivAssign, Mul, MulAssign};

/// Denominator representation that decomposes an integer as a product of the
/// first `NUM_PRIMES` primes (up to 2^16), multiplied by a regular big integer
/// when that's not sufficient.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DenomSparseU16<const NUM_PRIMES: usize> {
    // Invariant: the representation is in canonical form, i.e. powers of small primes must
    // saturate the primes array before overflowing into the remainder.
    primes: Vec<(u16, u16)>,
    remainder: Option<BigUint>,
}

/// Denominator representation that decomposes an integer as a product of the
/// first 6542 primes (up to 0xfff1), multiplied by a regular big integer when
/// that's not sufficient.
pub type DenomSparse6542 = DenomSparseU16<6542>;

impl<const NUM_PRIMES: usize> DenomRef<DenomSparseU16<NUM_PRIMES>> for &DenomSparseU16<NUM_PRIMES> {}

impl<const NUM_PRIMES: usize> Denom for DenomSparseU16<NUM_PRIMES> {
    const ONE: Self = Self {
        primes: Vec::new(),
        remainder: None,
    };

    fn into_biguint(self) -> BigUint {
        self.into()
    }

    fn to_biguint(&self) -> BigUint {
        self.into()
    }

    fn normalize(lnum: &mut BigInt, rnum: &mut BigInt, ldenom: &Self, rdenom: &Self) -> Self {
        let mut primes = Vec::new();
        let mut ltmp = 1_usize;
        let mut rtmp = 1_usize;

        for (p, (lcount, rcount)) in Zip(
            ldenom.primes.iter().copied().peekable(),
            rdenom.primes.iter().copied().peekable(),
        ) {
            match lcount.cmp(&rcount) {
                Ordering::Equal => {
                    primes.push((p, lcount));
                }
                Ordering::Less => {
                    Self::accum_pow(lnum, &mut ltmp, p, rcount - lcount);
                    primes.push((p, rcount));
                }
                Ordering::Greater => {
                    Self::accum_pow(rnum, &mut rtmp, p, lcount - rcount);
                    primes.push((p, lcount));
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

        'outer: for (p, i) in self.primes.iter_mut() {
            let p = BigInt::from(*p);
            while *i != 0 {
                let (quo, rem) = num.div_rem(&p);
                if !rem.is_zero() {
                    break;
                }
                *num = quo;
                *i -= 1;
                if num.magnitude().is_one() {
                    break 'outer;
                }
            }
        }
        self.primes.retain(|(_, i)| *i != 0);
    }
}

impl<const NUM_PRIMES: usize> DenomSparseU16<NUM_PRIMES> {
    fn decompose(mut x: BigUint) -> Self {
        const {
            assert!(NUM_PRIMES <= ODD_PRIMES.len() + 1);
        }

        let bits = x.bits();
        if bits <= 8 {
            return Self::decompose_u8(x.try_into().unwrap());
        } else if bits <= 16 {
            return Self::decompose_u16(x.try_into().unwrap());
        } else if bits <= 32 {
            return Self::decompose_u32(x.try_into().unwrap());
        } else if bits <= 64 {
            return Self::decompose_u64(x.try_into().unwrap());
        } else if bits <= 128 {
            return Self::decompose_u128(x.try_into().unwrap());
        }

        let mut primes = Vec::new();

        let mut count2 = x.trailing_zeros().unwrap();
        if count2 != 0 {
            x >>= count2;
            if count2 <= u16::MAX as u64 {
                primes.push((2, count2 as u16));
                count2 = 0;
            } else {
                primes.push((2, u16::MAX));
                count2 -= u16::MAX as u64;
            }
        }

        'outer: for &p in ODD_PRIMES.iter().take(NUM_PRIMES - 1) {
            let bigp = BigUint::from(p);
            let mut count = 0;
            while count != u16::MAX {
                let (quo, rem) = x.div_rem(&bigp);
                if !rem.is_zero() {
                    break;
                }
                x = quo;

                count += 1;
                if x.is_one() {
                    primes.push((p, count));
                    break 'outer;
                }
            }
            if count != 0 {
                primes.push((p, count));
            }
        }

        let remainder = if x.is_one() && count2 == 0 {
            None
        } else {
            Some(x << count2)
        };
        Self { primes, remainder }
    }

    fn decompose_known_factors(
        mut primes: Vec<(u16, u16)>,
        p: u16,
        mut pcount: u16,
        factor_indices: impl Iterator<Item = (u16, u16)>,
    ) -> Self {
        let mut remainder = 1_usize;
        for (index, mut qcount) in factor_indices {
            let index = index as usize;
            let q = ODD_PRIMES[index - 1];

            if pcount != 0 {
                if p == q {
                    qcount += pcount;
                } else {
                    primes.push((p, pcount));
                }
                pcount = 0;
            }

            if index < NUM_PRIMES {
                primes.push((q, qcount));
            } else {
                let q = q as usize;
                for _ in 0..qcount {
                    remainder *= q;
                }
            }
        }

        let remainder = if remainder == 1 {
            None
        } else {
            Some(remainder.into())
        };
        Self { primes, remainder }
    }

    fn decompose_u8(mut x: u8) -> Self {
        let mut primes = Vec::new();

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes.push((2, count2 as u16));
        }

        // 8-bit integers always fit in the look-up table.
        Self::decompose_known_factors(
            primes,
            0,
            0,
            known_odd_prime_factor_indices(x.into()).unwrap(),
        )
    }

    fn decompose_u16(mut x: u16) -> Self {
        const {
            assert!(NUM_PRIMES <= ODD_PRIMES.len() + 1);
        }

        let mut primes = Vec::new();

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes.push((2, count2 as u16));
        }

        'outer: for (i, &p) in ODD_PRIMES.iter().enumerate().take(NUM_PRIMES - 1) {
            let divider = OddDivider {
                divisor: p,
                multiplier: ODD_PRIME_DIVIDERS_U16[i],
                shift: p.ilog2(),
            };
            let mut count = 0;
            while x != 1 {
                // Use look-up table as soon as the remainder is small enough.
                if let Some(iter) = known_odd_prime_factor_indices(x.into()) {
                    return Self::decompose_known_factors(primes, p, count, iter);
                }

                let (quo, rem) = divider.div_rem(x);
                if rem != 0 {
                    break;
                }
                x = quo;

                count += 1;
                if x == 1 {
                    primes.push((p, count));
                    break 'outer;
                }
            }
            if count != 0 {
                primes.push((p, count));
            }
        }

        let remainder = if x == 1 { None } else { Some(x.into()) };
        Self { primes, remainder }
    }

    fn decompose_u32(mut x: u32) -> Self {
        const {
            assert!(NUM_PRIMES <= ODD_PRIMES.len() + 1);
        }

        let mut primes = Vec::new();

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes.push((2, count2 as u16));
        }

        'outer: for (i, &p) in ODD_PRIMES.iter().enumerate().take(NUM_PRIMES - 1) {
            let bigp = p as u32;
            let divider = OddDivider {
                divisor: bigp,
                multiplier: ODD_PRIME_DIVIDERS_U32[i],
                shift: p.ilog2(),
            };
            let mut count = 0;
            while x != 1 {
                // Use look-up table as soon as the remainder is small enough.
                if let Ok(xx) = x.try_into()
                    && let Some(iter) = known_odd_prime_factor_indices(xx)
                {
                    return Self::decompose_known_factors(primes, p, count, iter);
                }

                let (quo, rem) = divider.div_rem(x);
                if rem != 0 {
                    break;
                }
                x = quo;

                count += 1;
                if x == 1 {
                    primes.push((p, count));
                    break 'outer;
                }
            }
            if count != 0 {
                primes.push((p, count));
            }
        }

        let remainder = if x == 1 { None } else { Some(x.into()) };
        Self { primes, remainder }
    }

    fn decompose_u64(mut x: u64) -> Self {
        const {
            assert!(NUM_PRIMES <= ODD_PRIMES.len() + 1);
        }

        let mut primes = Vec::new();

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes.push((2, count2 as u16));
        }

        'outer: for (i, &p) in ODD_PRIMES.iter().enumerate().take(NUM_PRIMES - 1) {
            let bigp = p as u64;
            let divider = OddDivider {
                divisor: bigp,
                multiplier: ODD_PRIME_DIVIDERS_U64[i],
                shift: p.ilog2(),
            };
            let mut count = 0;
            while x != 1 {
                // Use look-up table as soon as the remainder is small enough.
                if let Ok(xx) = x.try_into()
                    && let Some(iter) = known_odd_prime_factor_indices(xx)
                {
                    return Self::decompose_known_factors(primes, p, count, iter);
                }

                let (quo, rem) = divider.div_rem(x);
                if rem != 0 {
                    break;
                }
                x = quo;

                count += 1;
                if x == 1 {
                    primes.push((p, count));
                    break 'outer;
                }
            }
            if count != 0 {
                primes.push((p, count));
            }
        }

        let remainder = if x == 1 { None } else { Some(x.into()) };
        Self { primes, remainder }
    }

    fn decompose_u128(mut x: u128) -> Self {
        const {
            assert!(NUM_PRIMES <= ODD_PRIMES.len() + 1);
        }

        let mut primes = Vec::new();

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes.push((2, count2 as u16));
        }

        'outer: for &p in ODD_PRIMES.iter().take(NUM_PRIMES - 1) {
            let bigp = p as u128;
            let mut count = 0;
            while x != 1 {
                // Use look-up table as soon as the remainder is small enough.
                if let Ok(xx) = x.try_into()
                    && let Some(iter) = known_odd_prime_factor_indices(xx)
                {
                    return Self::decompose_known_factors(primes, p, count, iter);
                }

                if !x.is_multiple_of(bigp) {
                    break;
                }
                x /= bigp;

                count += 1;
                if x == 1 {
                    primes.push((p, count));
                    break 'outer;
                }
            }
            if count != 0 {
                primes.push((p, count));
            }
        }

        let remainder = if x == 1 { None } else { Some(x.into()) };
        Self { primes, remainder }
    }

    fn decompose_usize(mut x: usize) -> Self {
        const {
            assert!(NUM_PRIMES <= ODD_PRIMES.len() + 1);
        }

        let mut primes = Vec::new();

        let count2 = x.trailing_zeros();
        if count2 != 0 {
            x >>= count2;
            primes.push((2, count2 as u16));
        }

        'outer: for &p in ODD_PRIMES.iter().take(NUM_PRIMES - 1) {
            let bigp = p as usize;
            let mut count = 0;
            while x != 1 {
                // Use look-up table as soon as the remainder is small enough.
                if let Some(iter) = known_odd_prime_factor_indices(x) {
                    return Self::decompose_known_factors(primes, p, count, iter);
                }

                if !x.is_multiple_of(bigp) {
                    break;
                }
                x /= bigp;

                count += 1;
                if x == 1 {
                    primes.push((p, count));
                    break 'outer;
                }
            }
            if count != 0 {
                primes.push((p, count));
            }
        }

        let remainder = if x == 1 { None } else { Some(x.into()) };
        Self { primes, remainder }
    }

    /// Decomposes the given remainder using only primes whose bit mask is set.
    fn decompose_mask(remainder: &mut Option<BigUint>, primes: &mut Vec<(u16, u16)>, mask: &[u16]) {
        let x: &mut BigUint = match remainder {
            None => return,
            Some(x) => x,
        };

        let mut new_primes = Vec::new();
        'outer: for (p, (mut count, flag)) in Zip(
            primes.iter().copied().peekable(),
            mask.iter().copied().map(|p| (p, true)).peekable(),
        ) {
            if flag && count != u16::MAX {
                let bigp = BigUint::from(p);
                while count != u16::MAX {
                    let (quo, rem) = x.div_rem(&bigp);
                    if !rem.is_zero() {
                        break;
                    }
                    *x = quo;

                    count += 1;
                    if x.is_one() {
                        new_primes.push((p, count));
                        break 'outer;
                    }
                }
            }
            new_primes.push((p, count));
        }
        *primes = new_primes;

        if x.is_one() {
            *remainder = None;
        }
    }

    /// Computes `prime.pow(exponent)` and multiplies it into the accumulated
    /// `(numerator, tmp)`.
    fn accum_pow(numerator: &mut BigInt, tmp: &mut usize, prime: u16, exponent: u16) {
        let prime = prime as usize;
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

    fn accum_pow_option(
        remainder: &mut Option<BigUint>,
        tmp: &mut usize,
        prime: u16,
        exponent: u16,
    ) {
        let prime = prime as usize;
        for _ in 0..exponent {
            match tmp.checked_mul(prime) {
                Some(prod) => *tmp = prod,
                None => {
                    match remainder {
                        None => *remainder = Some(BigUint::from(*tmp)),
                        Some(r) => *r *= *tmp,
                    };
                    *tmp = prime;
                }
            }
        }
    }

    fn pow_primes(primes: &mut [(u16, u16)], exponent: u32, remainder: &mut Option<BigUint>) {
        for (p, i) in primes.iter_mut() {
            let product = (*i as u32).strict_mul(exponent);
            if product <= u16::MAX as u32 {
                *i = product as u16;
            } else {
                *i = u16::MAX;
                let factor = BigUint::from(*p).pow(product - u16::MAX as u32);
                match remainder {
                    None => *remainder = Some(factor),
                    Some(r) => *r *= factor,
                };
            }
        }
    }

    fn mul_primes(
        lhs: &[(u16, u16)],
        rhs: &[(u16, u16)],
        remainder: &mut Option<BigUint>,
    ) -> Vec<(u16, u16)> {
        let mut primes = Vec::new();

        for (p, (lcount, rcount)) in Zip(
            lhs.iter().copied().peekable(),
            rhs.iter().copied().peekable(),
        ) {
            let sum = lcount as u32 + rcount as u32;
            if sum <= u16::MAX as u32 {
                primes.push((p, sum as u16));
            } else {
                primes.push((p, u16::MAX));
                let factor = BigUint::from(p).pow(sum - u16::MAX as u32);
                match remainder {
                    None => *remainder = Some(factor),
                    Some(r) => *r *= factor,
                };
            }
        }
        primes
    }

    fn div_primes(
        num: &[(u16, u16)],
        denom: &[(u16, u16)],
        mut remainder: Option<BigUint>,
    ) -> (Self, Vec<u16>) {
        let mut primes = Vec::new();
        let mut mask = Vec::new();
        let mut tmp = 1_usize;

        for (p, (num_count, denom_count)) in Zip(
            num.iter().copied().peekable(),
            denom.iter().copied().peekable(),
        ) {
            if num_count == u16::MAX {
                mask.push(p);
            }
            match num_count.cmp(&denom_count) {
                Ordering::Greater => {
                    primes.push((p, num_count - denom_count));
                }
                Ordering::Equal => (),
                Ordering::Less => {
                    Self::accum_pow_option(&mut remainder, &mut tmp, p, denom_count - num_count);
                    let factor = BigUint::from(p).pow(denom_count as u32 - num_count as u32);
                    remainder = match remainder {
                        None => Some(factor),
                        Some(r) => Some(r * factor),
                    };
                }
            }
        }

        if tmp != 1 {
            match remainder.as_mut() {
                None => remainder = Some(BigUint::from(tmp)),
                Some(r) => *r *= tmp,
            };
        }
        (Self { primes, remainder }, mask)
    }
}

impl<const NUM_PRIMES: usize> From<u8> for DenomSparseU16<NUM_PRIMES> {
    fn from(value: u8) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_u8(value)
    }
}

impl<const NUM_PRIMES: usize> From<u16> for DenomSparseU16<NUM_PRIMES> {
    fn from(value: u16) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_u16(value)
    }
}

impl<const NUM_PRIMES: usize> From<u32> for DenomSparseU16<NUM_PRIMES> {
    fn from(value: u32) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_u32(value)
    }
}

impl<const NUM_PRIMES: usize> From<u64> for DenomSparseU16<NUM_PRIMES> {
    fn from(value: u64) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_u64(value)
    }
}

impl<const NUM_PRIMES: usize> From<u128> for DenomSparseU16<NUM_PRIMES> {
    fn from(value: u128) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_u128(value)
    }
}

impl<const NUM_PRIMES: usize> From<usize> for DenomSparseU16<NUM_PRIMES> {
    fn from(value: usize) -> Self {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose_usize(value)
    }
}

impl<const NUM_PRIMES: usize> From<BigUint> for DenomSparseU16<NUM_PRIMES> {
    fn from(value: BigUint) -> Self {
        if value.is_zero() {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose(value)
    }
}

impl<const NUM_PRIMES: usize> From<&BigUint> for DenomSparseU16<NUM_PRIMES> {
    fn from(value: &BigUint) -> Self {
        if value.is_zero() {
            panic!("Attempted to create a denominator of zero");
        }
        Self::decompose(value.clone())
    }
}

impl<const NUM_PRIMES: usize> From<DenomSparseU16<NUM_PRIMES>> for BigUint {
    fn from(value: DenomSparseU16<NUM_PRIMES>) -> BigUint {
        let mut result = match value.remainder {
            Some(x) => x,
            None => BigUint::ONE,
        };
        let mut tmp = 1_usize;
        for (p, count) in value.primes.into_iter() {
            let p = p as usize;
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

impl<const NUM_PRIMES: usize> From<&DenomSparseU16<NUM_PRIMES>> for BigUint {
    fn from(value: &DenomSparseU16<NUM_PRIMES>) -> BigUint {
        let mut result = match &value.remainder {
            Some(x) => x.clone(),
            None => BigUint::ONE,
        };
        let mut tmp = 1_usize;
        for &(p, count) in value.primes.iter() {
            let p = p as usize;
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

impl<const NUM_PRIMES: usize> One for DenomSparseU16<NUM_PRIMES> {
    fn one() -> Self {
        Self::ONE
    }
}

impl<const NUM_PRIMES: usize> num_traits::Pow<u32> for DenomSparseU16<NUM_PRIMES> {
    type Output = Self;

    fn pow(mut self, rhs: u32) -> Self {
        if rhs == 0 {
            return Self::ONE;
        }
        Self::pow_primes(&mut self.primes, rhs, &mut self.remainder);
        self
    }
}

impl<const NUM_PRIMES: usize> num_traits::Pow<u32> for &DenomSparseU16<NUM_PRIMES> {
    type Output = DenomSparseU16<NUM_PRIMES>;

    fn pow(self, rhs: u32) -> DenomSparseU16<NUM_PRIMES> {
        if rhs == 0 {
            return DenomSparseU16::ONE;
        }
        let mut primes = self.primes.clone();
        let mut remainder = self.remainder.clone();
        DenomSparseU16::<NUM_PRIMES>::pow_primes(&mut primes, rhs, &mut remainder);
        DenomSparseU16 { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Mul for DenomSparseU16<NUM_PRIMES> {
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

impl<const NUM_PRIMES: usize> Mul<&Self> for DenomSparseU16<NUM_PRIMES> {
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

impl<const NUM_PRIMES: usize> Mul for &DenomSparseU16<NUM_PRIMES> {
    type Output = DenomSparseU16<NUM_PRIMES>;

    fn mul(self, rhs: Self) -> DenomSparseU16<NUM_PRIMES> {
        let mut remainder = match (&self.remainder, &rhs.remainder) {
            (None, None) => None,
            (None, Some(r)) => Some(r.clone()),
            (Some(l), None) => Some(l.clone()),
            (Some(l), Some(r)) => Some(l * r),
        };
        let primes =
            DenomSparseU16::<NUM_PRIMES>::mul_primes(&self.primes, &rhs.primes, &mut remainder);
        DenomSparseU16 { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Mul<DenomSparseU16<NUM_PRIMES>> for &DenomSparseU16<NUM_PRIMES> {
    type Output = DenomSparseU16<NUM_PRIMES>;

    fn mul(self, rhs: DenomSparseU16<NUM_PRIMES>) -> DenomSparseU16<NUM_PRIMES> {
        let mut remainder = match (&self.remainder, rhs.remainder) {
            (None, None) => None,
            (None, Some(r)) => Some(r),
            (Some(l), None) => Some(l.clone()),
            (Some(l), Some(r)) => Some(l * r),
        };
        let primes =
            DenomSparseU16::<NUM_PRIMES>::mul_primes(&self.primes, &rhs.primes, &mut remainder);
        DenomSparseU16 { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> MulAssign for DenomSparseU16<NUM_PRIMES> {
    fn mul_assign(&mut self, rhs: Self) {
        match (&mut self.remainder, rhs.remainder) {
            (_, None) => (),
            (None, Some(r)) => self.remainder = Some(r),
            (Some(l), Some(r)) => *l *= r,
        };
        self.primes = Self::mul_primes(&self.primes, &rhs.primes, &mut self.remainder);
    }
}

impl<const NUM_PRIMES: usize> MulAssign<&Self> for DenomSparseU16<NUM_PRIMES> {
    fn mul_assign(&mut self, rhs: &Self) {
        match (&mut self.remainder, &rhs.remainder) {
            (_, None) => (),
            (None, Some(r)) => self.remainder = Some(r.clone()),
            (Some(l), Some(r)) => *l *= r,
        };
        self.primes = Self::mul_primes(&self.primes, &rhs.primes, &mut self.remainder);
    }
}

impl<const NUM_PRIMES: usize> Div for DenomSparseU16<NUM_PRIMES> {
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

        Self::decompose_mask(&mut remainder, &mut primes, &mask);
        Self { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Div<&Self> for DenomSparseU16<NUM_PRIMES> {
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

        Self::decompose_mask(&mut remainder, &mut primes, &mask);
        Self { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Div for &DenomSparseU16<NUM_PRIMES> {
    type Output = DenomSparseU16<NUM_PRIMES>;

    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div(self, rhs: Self) -> DenomSparseU16<NUM_PRIMES> {
        let (
            DenomSparseU16 {
                mut primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = DenomSparseU16::<NUM_PRIMES>::div_primes(
            &self.primes,
            &rhs.primes,
            rhs.remainder.clone(),
        );

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

        DenomSparseU16::<NUM_PRIMES>::decompose_mask(&mut remainder, &mut primes, &mask);
        DenomSparseU16 { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> Div<DenomSparseU16<NUM_PRIMES>> for &DenomSparseU16<NUM_PRIMES> {
    type Output = DenomSparseU16<NUM_PRIMES>;

    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div(self, rhs: DenomSparseU16<NUM_PRIMES>) -> DenomSparseU16<NUM_PRIMES> {
        let (
            DenomSparseU16 {
                mut primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = DenomSparseU16::<NUM_PRIMES>::div_primes(&self.primes, &rhs.primes, rhs.remainder);

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

        DenomSparseU16::<NUM_PRIMES>::decompose_mask(&mut remainder, &mut primes, &mask);
        DenomSparseU16 { primes, remainder }
    }
}

impl<const NUM_PRIMES: usize> DivAssign for DenomSparseU16<NUM_PRIMES> {
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

        Self::decompose_mask(&mut self.remainder, &mut self.primes, &mask);
    }
}

impl<const NUM_PRIMES: usize> DivAssign<&Self> for DenomSparseU16<NUM_PRIMES> {
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

        Self::decompose_mask(&mut self.remainder, &mut self.primes, &mask);
    }
}

impl<const NUM_PRIMES: usize> Mul<DenomSparseU16<NUM_PRIMES>> for BigInt {
    type Output = Self;

    fn mul(mut self, rhs: DenomSparseU16<NUM_PRIMES>) -> Self {
        let mut tmp = 1_usize;
        for (p, count) in rhs.primes.into_iter() {
            DenomSparseU16::<NUM_PRIMES>::accum_pow(&mut self, &mut tmp, p, count);
        }

        self *= tmp;
        if let Some(remainder) = rhs.remainder {
            self *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        self
    }
}

impl<const NUM_PRIMES: usize> Mul<&DenomSparseU16<NUM_PRIMES>> for BigInt {
    type Output = Self;

    fn mul(mut self, rhs: &DenomSparseU16<NUM_PRIMES>) -> Self {
        let mut tmp = 1_usize;
        for &(p, count) in rhs.primes.iter() {
            DenomSparseU16::<NUM_PRIMES>::accum_pow(&mut self, &mut tmp, p, count);
        }

        self *= tmp;
        if let Some(remainder) = &rhs.remainder {
            self *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        self
    }
}

impl<const NUM_PRIMES: usize> Mul<DenomSparseU16<NUM_PRIMES>> for &BigInt {
    type Output = BigInt;

    fn mul(self, rhs: DenomSparseU16<NUM_PRIMES>) -> BigInt {
        let mut this = self.clone();
        let mut tmp = 1_usize;
        for (p, count) in rhs.primes.into_iter() {
            DenomSparseU16::<NUM_PRIMES>::accum_pow(&mut this, &mut tmp, p, count);
        }

        this *= tmp;
        if let Some(remainder) = rhs.remainder {
            this *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        this
    }
}

impl<const NUM_PRIMES: usize> Mul<&DenomSparseU16<NUM_PRIMES>> for &BigInt {
    type Output = BigInt;

    fn mul(self, rhs: &DenomSparseU16<NUM_PRIMES>) -> BigInt {
        let mut this = self.clone();
        let mut tmp = 1_usize;
        for &(p, count) in rhs.primes.iter() {
            DenomSparseU16::<NUM_PRIMES>::accum_pow(&mut this, &mut tmp, p, count);
        }

        this *= tmp;
        if let Some(remainder) = &rhs.remainder {
            this *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        this
    }
}

impl<const NUM_PRIMES: usize> Mul<BigInt> for DenomSparseU16<NUM_PRIMES> {
    type Output = BigInt;

    fn mul(self, mut rhs: BigInt) -> BigInt {
        let mut tmp = 1_usize;
        for (p, count) in self.primes.into_iter() {
            Self::accum_pow(&mut rhs, &mut tmp, p, count);
        }

        rhs *= tmp;
        if let Some(remainder) = self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        rhs
    }
}

impl<const NUM_PRIMES: usize> Mul<&BigInt> for DenomSparseU16<NUM_PRIMES> {
    type Output = BigInt;

    fn mul(self, rhs: &BigInt) -> BigInt {
        let mut rhs = rhs.clone();
        let mut tmp = 1_usize;
        for (p, count) in self.primes.into_iter() {
            Self::accum_pow(&mut rhs, &mut tmp, p, count);
        }

        rhs *= tmp;
        if let Some(remainder) = self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        rhs
    }
}

impl<const NUM_PRIMES: usize> Mul<BigInt> for &DenomSparseU16<NUM_PRIMES> {
    type Output = BigInt;

    fn mul(self, mut rhs: BigInt) -> BigInt {
        let mut tmp = 1_usize;
        for &(p, count) in self.primes.iter() {
            DenomSparseU16::<NUM_PRIMES>::accum_pow(&mut rhs, &mut tmp, p, count);
        }

        rhs *= tmp;
        if let Some(remainder) = &self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        rhs
    }
}

impl<const NUM_PRIMES: usize> Mul<&BigInt> for &DenomSparseU16<NUM_PRIMES> {
    type Output = BigInt;

    fn mul(self, rhs: &BigInt) -> BigInt {
        let mut rhs = rhs.clone();
        let mut tmp = 1_usize;
        for &(p, count) in self.primes.iter() {
            DenomSparseU16::<NUM_PRIMES>::accum_pow(&mut rhs, &mut tmp, p, count);
        }

        rhs *= tmp;
        if let Some(remainder) = &self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        rhs
    }
}

impl<const NUM_PRIMES: usize> MulAssign<DenomSparseU16<NUM_PRIMES>> for BigInt {
    fn mul_assign(&mut self, rhs: DenomSparseU16<NUM_PRIMES>) {
        let mut tmp = 1_usize;
        for (p, count) in rhs.primes.into_iter() {
            DenomSparseU16::<NUM_PRIMES>::accum_pow(self, &mut tmp, p, count);
        }

        *self *= tmp;
        if let Some(remainder) = rhs.remainder {
            *self *= BigInt::from_biguint(Sign::Plus, remainder);
        }
    }
}

impl<const NUM_PRIMES: usize> MulAssign<&DenomSparseU16<NUM_PRIMES>> for BigInt {
    fn mul_assign(&mut self, rhs: &DenomSparseU16<NUM_PRIMES>) {
        let mut tmp = 1_usize;
        for &(p, count) in rhs.primes.iter() {
            DenomSparseU16::<NUM_PRIMES>::accum_pow(self, &mut tmp, p, count);
        }

        *self *= tmp;
        if let Some(remainder) = &rhs.remainder {
            *self *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
    }
}

struct Zip<A, B>(A, B);

impl<U, V, A, B> Iterator for Zip<Peekable<A>, Peekable<B>>
where
    A: Iterator<Item = (u16, U)>,
    B: Iterator<Item = (u16, V)>,
    U: Default,
    V: Default,
{
    type Item = (u16, (U, V));

    fn next(&mut self) -> Option<Self::Item> {
        match (self.0.peek(), self.1.peek()) {
            (None, None) => None,
            (None, Some(_)) => {
                let (p, v) = self.1.next().unwrap();
                Some((p, (Default::default(), v)))
            }
            (Some(_), None) => {
                let (p, u) = self.0.next().unwrap();
                Some((p, (u, Default::default())))
            }
            (Some((p, _)), Some((q, _))) => match p.cmp(q) {
                Ordering::Less => {
                    let (p, u) = self.0.next().unwrap();
                    Some((p, (u, Default::default())))
                }
                Ordering::Greater => {
                    let (q, v) = self.1.next().unwrap();
                    Some((q, (Default::default(), v)))
                }
                Ordering::Equal => {
                    let (p, u) = self.0.next().unwrap();
                    let (_, v) = self.1.next().unwrap();
                    Some((p, (u, v)))
                }
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! tests {
        (
            $mod:ident,
            $num_primes:expr,
            $( $case:ident ,)*
        ) => {
            mod $mod {
                $(
                    #[test]
                    fn $case() {
                        super::$case::<$num_primes>();
                    }
                )*
            }
        };
    }

    macro_rules! all_tests {
        (
            $mod:ident,
            $num_primes:expr
        ) => {
            tests!(
                $mod,
                $num_primes,
                test_decompose_to_biguint,
                test_decompose_small_prime_power,
                test_mul_prime_powers,
                test_div_prime_powers,
                test_product,
                test_normalize,
                test_gcd_reduce,
            );
        };
    }

    all_tests!(denom24, 24);
    all_tests!(denom6542, 6542);

    fn test_decompose_to_biguint<const NUM_PRIMES: usize>() {
        for i in 1_usize..=(1 << 20) {
            let bigi = BigUint::from(i);
            let x = DenomSparseU16::<NUM_PRIMES>::from(&bigi);
            assert_eq!(x.to_biguint(), bigi);
        }
    }

    #[test]
    fn test_decompose_no_remainder() {
        for i in 1_usize..=65536 {
            let bigi = BigUint::from(i);
            let x = DenomSparse6542::from(&bigi);
            assert_eq!(x.remainder, None, "Remainder for {i}");
        }
    }

    #[test]
    fn test_decompose_known_values() {
        assert_eq!(
            DenomSparse6542::from(0xfff1_u16),
            DenomSparse6542 {
                primes: vec![(0xfff1, 1)],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(0xffff_fffb_u32),
            DenomSparse6542 {
                primes: vec![],
                remainder: Some(0xffff_fffb_u32.into()),
            }
        );
        assert_eq!(
            DenomSparse6542::from(0xffff_ffff_ffff_ffc5_u64),
            DenomSparse6542 {
                primes: vec![],
                remainder: Some(0xffff_ffff_ffff_ffc5_u64.into()),
            }
        );
        assert_eq!(
            DenomSparse6542::from(0xffff_ffff_ffff_ffff_ffff_ffff_ffff_ff61_u128),
            DenomSparse6542 {
                primes: vec![],
                remainder: Some(0xffff_ffff_ffff_ffff_ffff_ffff_ffff_ff61_u128.into()),
            }
        );
        assert_eq!(
            DenomSparse6542::from(BigUint::from_slice(&[
                0xffff_ff43,
                0xffff_ffff,
                0xffff_ffff,
                0xffff_ffff,
                0xffff_ffff,
                0xffff_ffff,
                0xffff_ffff,
                0xffff_ffff,
            ])),
            DenomSparse6542 {
                primes: vec![],
                remainder: Some(BigUint::from_slice(&[
                    0xffff_ff43,
                    0xffff_ffff,
                    0xffff_ffff,
                    0xffff_ffff,
                    0xffff_ffff,
                    0xffff_ffff,
                    0xffff_ffff,
                    0xffff_ffff,
                ])),
            }
        );

        assert_eq!(
            DenomSparse6542::from(2u8 * 3 * 5 * 7),
            DenomSparse6542 {
                primes: vec![(2, 1), (3, 1), (5, 1), (7, 1)],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(2u16 * 3 * 5 * 7 * 11 * 13),
            DenomSparse6542 {
                primes: vec![(2, 1), (3, 1), (5, 1), (7, 1), (11, 1), (13, 1)],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(2u32 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23),
            DenomSparse6542 {
                primes: vec![
                    (2, 1),
                    (3, 1),
                    (5, 1),
                    (7, 1),
                    (11, 1),
                    (13, 1),
                    (17, 1),
                    (19, 1),
                    (23, 1)
                ],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(
                2u64 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 31 * 37 * 41 * 43 * 47 * 53
            ),
            DenomSparse6542 {
                primes: vec![
                    (2, 1),
                    (3, 1),
                    (5, 1),
                    (7, 1),
                    (11, 1),
                    (13, 1),
                    (17, 1),
                    (19, 1),
                    (23, 1),
                    (31, 1),
                    (37, 1),
                    (41, 1),
                    (43, 1),
                    (47, 1),
                    (53, 1)
                ],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(
                2u128
                    * 3
                    * 5
                    * 7
                    * 11
                    * 13
                    * 17
                    * 19
                    * 23
                    * 31
                    * 37
                    * 41
                    * 43
                    * 47
                    * 53
                    * 59
                    * 61
                    * 67
                    * 71
                    * 73
                    * 79
                    * 83
                    * 89
                    * 97
                    * 101,
            ),
            DenomSparse6542 {
                primes: vec![
                    (2, 1),
                    (3, 1),
                    (5, 1),
                    (7, 1),
                    (11, 1),
                    (13, 1),
                    (17, 1),
                    (19, 1),
                    (23, 1),
                    (31, 1),
                    (37, 1),
                    (41, 1),
                    (43, 1),
                    (47, 1),
                    (53, 1),
                    (59, 1),
                    (61, 1),
                    (67, 1),
                    (71, 1),
                    (73, 1),
                    (79, 1),
                    (83, 1),
                    (89, 1),
                    (97, 1),
                    (101, 1),
                ],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(
                BigUint::try_from(
                    BigInt::from(2)
                        * 3
                        * 5
                        * 7
                        * 11
                        * 13
                        * 17
                        * 19
                        * 23
                        * 31
                        * 37
                        * 41
                        * 43
                        * 47
                        * 53
                        * 59
                        * 61
                        * 67
                        * 71
                        * 73
                        * 79
                        * 83
                        * 89
                        * 97
                        * 101
                        * 103
                        * 107
                        * 109
                        * 113
                        * 127
                        * 131
                        * 137
                        * 139
                        * 149
                        * 151
                        * 157
                        * 163
                        * 167
                        * 173
                        * 179
                        * 181
                        * 191
                        * 193
                        * 197
                        * 199
                        * 211
                        * 223
                        * 227
                        * 229
                        * 233
                )
                .unwrap(),
            ),
            DenomSparse6542 {
                primes: vec![
                    (2, 1),
                    (3, 1),
                    (5, 1),
                    (7, 1),
                    (11, 1),
                    (13, 1),
                    (17, 1),
                    (19, 1),
                    (23, 1),
                    (31, 1),
                    (37, 1),
                    (41, 1),
                    (43, 1),
                    (47, 1),
                    (53, 1),
                    (59, 1),
                    (61, 1),
                    (67, 1),
                    (71, 1),
                    (73, 1),
                    (79, 1),
                    (83, 1),
                    (89, 1),
                    (97, 1),
                    (101, 1),
                    (103, 1),
                    (107, 1),
                    (109, 1),
                    (113, 1),
                    (127, 1),
                    (131, 1),
                    (137, 1),
                    (139, 1),
                    (149, 1),
                    (151, 1),
                    (157, 1),
                    (163, 1),
                    (167, 1),
                    (173, 1),
                    (179, 1),
                    (181, 1),
                    (191, 1),
                    (193, 1),
                    (197, 1),
                    (199, 1),
                    (211, 1),
                    (223, 1),
                    (227, 1),
                    (229, 1),
                    (233, 1),
                ],
                remainder: None,
            }
        );

        assert_eq!(
            DenomSparse6542::from(0xfb_u16 * 0xf1),
            DenomSparse6542 {
                primes: vec![(0xf1, 1), (0xfb, 1)],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(0xfff1_u32 * 0xffef),
            DenomSparse6542 {
                primes: vec![(0xffef, 1), (0xfff1, 1)],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(0xfff1_u64 * 0xffef * 0xffd9 * 0xffc7),
            DenomSparse6542 {
                primes: vec![(0xffc7, 1), (0xffd9, 1), (0xffef, 1), (0xfff1, 1)],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(
                0xfff1_u128 * 0xffef * 0xffd9 * 0xffc7 * 0xffa9 * 0xffa7 * 0xff9d * 0xff8f,
            ),
            DenomSparse6542 {
                primes: vec![
                    (0xff8f, 1),
                    (0xff9d, 1),
                    (0xffa7, 1),
                    (0xffa9, 1),
                    (0xffc7, 1),
                    (0xffd9, 1),
                    (0xffef, 1),
                    (0xfff1, 1)
                ],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(
                BigUint::from(0xfff1_u16)
                    * 0xffef_u16
                    * 0xffd9_u16
                    * 0xffc7_u16
                    * 0xffa9_u16
                    * 0xffa7_u16
                    * 0xff9d_u16
                    * 0xff8f_u16
                    * 0xff8b_u16
                    * 0xff85_u16
                    * 0xff7f_u16
                    * 0xff71_u16
                    * 0xff65_u16
                    * 0xff5b_u16
                    * 0xff4d_u16
                    * 0xff49_u16,
            ),
            DenomSparse6542 {
                primes: vec![
                    (0xff49, 1),
                    (0xff4d, 1),
                    (0xff5b, 1),
                    (0xff65, 1),
                    (0xff71, 1),
                    (0xff7f, 1),
                    (0xff85, 1),
                    (0xff8b, 1),
                    (0xff8f, 1),
                    (0xff9d, 1),
                    (0xffa7, 1),
                    (0xffa9, 1),
                    (0xffc7, 1),
                    (0xffd9, 1),
                    (0xffef, 1),
                    (0xfff1, 1)
                ],
                remainder: None,
            }
        );

        assert_eq!(
            DenomSparse6542::from(BigUint::from(128_usize)),
            DenomSparse6542 {
                primes: vec![(2, 7)],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(BigUint::from(89_usize)),
            DenomSparse6542 {
                primes: vec![(89, 1)],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(BigUint::from(97_usize)),
            DenomSparse6542 {
                primes: vec![(97, 1)],
                remainder: None,
            }
        );
        assert_eq!(
            DenomSparse6542::from(BigUint::from(97000_usize)),
            DenomSparse6542 {
                primes: vec![(2, 3), (5, 3), (97, 1)],
                remainder: None,
            }
        );
    }

    fn test_decompose_small_prime_power<const NUM_PRIMES: usize>() {
        for p in std::iter::once(2).chain(ODD_PRIMES).take(NUM_PRIMES) {
            let bigp = BigUint::from(p);
            assert_eq!(
                DenomSparseU16::<NUM_PRIMES>::from(bigp.pow(100)),
                DenomSparseU16 {
                    primes: vec![(p, 100)],
                    remainder: None,
                }
            );
        }
    }

    fn test_mul_prime_powers<const NUM_PRIMES: usize>() {
        let p = BigUint::from(2u32);
        for a in 1..=256 {
            let denom_a = DenomSparseU16::<NUM_PRIMES>::from(p.pow(a));
            for b in 1..=256 {
                let denom_b = DenomSparseU16::<NUM_PRIMES>::from(p.pow(b));
                let denom_ab = DenomSparseU16::<NUM_PRIMES>::from(p.pow(a + b));
                assert_eq!(&denom_a * denom_b, denom_ab);
            }
        }
    }

    fn test_div_prime_powers<const NUM_PRIMES: usize>() {
        let p = BigUint::from(2u32);
        for a in 1..=256 {
            let denom_a = DenomSparseU16::<NUM_PRIMES>::from(p.pow(a));
            for b in 1..=256 {
                let denom_b = DenomSparseU16::<NUM_PRIMES>::from(p.pow(b));
                let denom_ab = DenomSparseU16::<NUM_PRIMES>::from(p.pow(a + b));
                assert_eq!(denom_ab / denom_b, denom_a);
            }
        }
    }

    fn test_product<const NUM_PRIMES: usize>() {
        let values = (100..200)
            .map(|i: usize| DenomSparseU16::<NUM_PRIMES>::from(BigUint::from(i)))
            .collect::<Vec<_>>();
        for (i, x) in values.iter().enumerate().map(|(i, x)| (i + 100, x)) {
            for (j, y) in values.iter().enumerate().map(|(j, y)| (j + 100, y)) {
                let z = x * y;
                assert_eq!(z, DenomSparseU16::<NUM_PRIMES>::from(BigUint::from(i * j)));

                for (p, (zcount, (xcount, ycount))) in Zip(
                    z.primes.into_iter().peekable(),
                    Zip(
                        x.primes.iter().copied().peekable(),
                        y.primes.iter().copied().peekable(),
                    )
                    .peekable(),
                ) {
                    assert_eq!(
                        zcount,
                        xcount + ycount,
                        "{zcount} != {xcount} + {ycount} for {p}"
                    );
                }
            }
        }
    }

    fn test_normalize<const NUM_PRIMES: usize>() {
        let values = (100..200)
            .map(|i: usize| DenomSparseU16::<NUM_PRIMES>::from(BigUint::from(i)))
            .collect::<Vec<_>>();
        for x in &values {
            for y in &values {
                let mut xnum = BigInt::one();
                let mut ynum = BigInt::one();
                let lcm = DenomSparseU16::<NUM_PRIMES>::normalize(&mut xnum, &mut ynum, x, y);
                let lcm_bigint = lcm.to_biguint();
                let xnum = xnum.to_biguint().unwrap();
                let ynum = ynum.to_biguint().unwrap();

                assert_eq!(xnum * x.to_biguint(), lcm_bigint);
                assert_eq!(ynum * y.to_biguint(), lcm_bigint);

                for (p, (lcm_count, (xcount, ycount))) in Zip(
                    lcm.primes.into_iter().peekable(),
                    Zip(
                        x.primes.iter().copied().peekable(),
                        y.primes.iter().copied().peekable(),
                    )
                    .peekable(),
                ) {
                    assert_eq!(
                        lcm_count,
                        std::cmp::max(xcount, ycount),
                        "{lcm_count} != max({xcount}, {ycount}) for {p}"
                    );
                }
            }
        }
    }

    fn test_gcd_reduce<const NUM_PRIMES: usize>() {
        let mut num = BigInt::from(-3 * 97);
        let mut denom = DenomSparseU16::<NUM_PRIMES>::from(BigUint::from(3u32 * 5 * 97));
        denom.gcd_reduce(&mut num);
        assert_eq!(num, BigInt::from(-1));
        assert_eq!(
            denom,
            DenomSparseU16::<NUM_PRIMES>::from(BigUint::from(5u32))
        );
    }
}
