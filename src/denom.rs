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

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};
use std::cmp::Ordering;
use std::ops::{Div, DivAssign, Mul, MulAssign};

/// Helper struct to represent the positive denominator of a rational number.
///
/// Under the hood, this decomposes the underlying integer as a product of small
/// primes, multiplied by a regular big integer when that's not sufficient.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Denom {
    // Invariant: the representation is in canonical form, i.e. powers of small primes must
    // saturate the primes array before overflowing into the remainder.
    primes: [u8; Self::NUM_PRIMES],
    remainder: Option<BigUint>,
}

impl Denom {
    const NUM_PRIMES: usize = 24;
    const PRIMES: [usize; Self::NUM_PRIMES] = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
    ];

    /// Constant value of 1.
    pub const ONE: Self = Denom {
        primes: [0; Self::NUM_PRIMES],
        remainder: None,
    };

    /// Converts this denominator into a big integer.
    pub fn into_biguint(self) -> BigUint {
        self.into()
    }

    /// Converts this denominator into a big integer.
    pub fn to_biguint(&self) -> BigUint {
        self.into()
    }

    fn decompose_now(mut x: BigUint) -> Self {
        let mut primes = [0; Self::NUM_PRIMES];
        'outer: for (i, &p) in Self::PRIMES.iter().enumerate() {
            let p = BigUint::from(p);
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

        let remainder = if x.is_one() { None } else { Some(x) };
        Self { primes, remainder }
    }

    fn decompose_u8(mut x: u8) -> Self {
        let mut primes = [0; Self::NUM_PRIMES];
        'outer: for (i, &p) in Self::PRIMES.iter().enumerate() {
            let p = p as u8;
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
        let mut primes = [0; Self::NUM_PRIMES];
        'outer: for (i, &p) in Self::PRIMES.iter().enumerate() {
            let p = p as u16;
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
        let mut primes = [0; Self::NUM_PRIMES];
        'outer: for (i, &p) in Self::PRIMES.iter().enumerate() {
            let p = p as u32;
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
        let mut primes = [0; Self::NUM_PRIMES];
        'outer: for (i, &p) in Self::PRIMES.iter().enumerate() {
            let p = p as u64;
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
        let mut primes = [0; Self::NUM_PRIMES];
        'outer: for (i, &p) in Self::PRIMES.iter().enumerate() {
            let p = p as u128;
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
        let mut primes = [0; Self::NUM_PRIMES];
        'outer: for (i, &p) in Self::PRIMES.iter().enumerate() {
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
        primes: &mut [u8; Self::NUM_PRIMES],
        mask: [bool; Self::NUM_PRIMES],
    ) {
        let x: &mut BigUint = match remainder {
            None => return,
            Some(x) => x,
        };

        'outer: for (i, &p) in Self::PRIMES.iter().enumerate() {
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

    /// Decomposes the given value into a denominator consisting of small prime
    /// factors.
    ///
    /// # Panics
    ///
    /// This function panics if the input isn't a product of small prime
    /// factors, up to the 24th prime number i.e. 89.
    pub const fn decompose_small(mut x: usize) -> Self {
        let mut primes = [0u8; Self::NUM_PRIMES];
        // TODO: use a for loop once supported in `const fn` context.
        let mut i = 0;
        while x > 1 && i < Self::NUM_PRIMES {
            let p = Self::PRIMES[i];
            while x.is_multiple_of(p) {
                x /= p;
                primes[i] = primes[i].checked_add(1).unwrap();
            }
            i += 1;
        }

        if x != 1 {
            panic!("Failed to decompose small integer into small prime factors.");
        }
        Self {
            primes,
            remainder: None,
        }
    }

    // TODO: Use std::array::from_fn when available in const contexts.
    const DECOMPOSED: [Self; 90] = [
        Self::decompose_small(1),
        Self::decompose_small(2),
        Self::decompose_small(3),
        Self::decompose_small(4),
        Self::decompose_small(5),
        Self::decompose_small(6),
        Self::decompose_small(7),
        Self::decompose_small(8),
        Self::decompose_small(9),
        Self::decompose_small(10),
        Self::decompose_small(11),
        Self::decompose_small(12),
        Self::decompose_small(13),
        Self::decompose_small(14),
        Self::decompose_small(15),
        Self::decompose_small(16),
        Self::decompose_small(17),
        Self::decompose_small(18),
        Self::decompose_small(19),
        Self::decompose_small(20),
        Self::decompose_small(21),
        Self::decompose_small(22),
        Self::decompose_small(23),
        Self::decompose_small(24),
        Self::decompose_small(25),
        Self::decompose_small(26),
        Self::decompose_small(27),
        Self::decompose_small(28),
        Self::decompose_small(29),
        Self::decompose_small(30),
        Self::decompose_small(31),
        Self::decompose_small(32),
        Self::decompose_small(33),
        Self::decompose_small(34),
        Self::decompose_small(35),
        Self::decompose_small(36),
        Self::decompose_small(37),
        Self::decompose_small(38),
        Self::decompose_small(39),
        Self::decompose_small(40),
        Self::decompose_small(41),
        Self::decompose_small(42),
        Self::decompose_small(43),
        Self::decompose_small(44),
        Self::decompose_small(45),
        Self::decompose_small(46),
        Self::decompose_small(47),
        Self::decompose_small(48),
        Self::decompose_small(49),
        Self::decompose_small(50),
        Self::decompose_small(51),
        Self::decompose_small(52),
        Self::decompose_small(53),
        Self::decompose_small(54),
        Self::decompose_small(55),
        Self::decompose_small(56),
        Self::decompose_small(57),
        Self::decompose_small(58),
        Self::decompose_small(59),
        Self::decompose_small(60),
        Self::decompose_small(61),
        Self::decompose_small(62),
        Self::decompose_small(63),
        Self::decompose_small(64),
        Self::decompose_small(65),
        Self::decompose_small(66),
        Self::decompose_small(67),
        Self::decompose_small(68),
        Self::decompose_small(69),
        Self::decompose_small(70),
        Self::decompose_small(71),
        Self::decompose_small(72),
        Self::decompose_small(73),
        Self::decompose_small(74),
        Self::decompose_small(75),
        Self::decompose_small(76),
        Self::decompose_small(77),
        Self::decompose_small(78),
        Self::decompose_small(79),
        Self::decompose_small(80),
        Self::decompose_small(81),
        Self::decompose_small(82),
        Self::decompose_small(83),
        Self::decompose_small(84),
        Self::decompose_small(85),
        Self::decompose_small(86),
        Self::decompose_small(87),
        Self::decompose_small(88),
        Self::decompose_small(89),
        Self::decompose_small(90),
    ];

    /// Returns the least common multiple of two denominators, adjusting the
    /// numerators accordingly.
    pub fn normalize(lnum: &mut BigInt, rnum: &mut BigInt, ldenom: &Self, rdenom: &Self) -> Self {
        let mut primes = [0; Self::NUM_PRIMES];
        let mut ltmp = 1_usize;
        let mut rtmp = 1_usize;
        for (i, &p) in Self::PRIMES.iter().enumerate() {
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
        Denom { primes, remainder }
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
        this: &[u8; Self::NUM_PRIMES],
        exponent: u32,
        remainder: &mut Option<BigUint>,
    ) -> [u8; Self::NUM_PRIMES] {
        let mut primes = [0; Self::NUM_PRIMES];
        for (i, &p) in Self::PRIMES.iter().enumerate() {
            let product = this[i] as u32 * exponent;
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
        lhs: &[u8; Self::NUM_PRIMES],
        rhs: &[u8; Self::NUM_PRIMES],
        remainder: &mut Option<BigUint>,
    ) -> [u8; Self::NUM_PRIMES] {
        let mut primes = [0; Self::NUM_PRIMES];
        for (i, &p) in Self::PRIMES.iter().enumerate() {
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
        num: &[u8; Self::NUM_PRIMES],
        denom: &[u8; Self::NUM_PRIMES],
        mut remainder: Option<BigUint>,
    ) -> (Self, [bool; Self::NUM_PRIMES]) {
        let mut primes = [0; Self::NUM_PRIMES];
        let mut mask = [false; Self::NUM_PRIMES];
        for (i, &p) in Self::PRIMES.iter().enumerate() {
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
        (Denom { primes, remainder }, mask)
    }

    /// Reduces this denominator together with the given numerator so that their
    /// GCD is one.
    pub fn gcd_reduce(&mut self, num: &mut BigInt) {
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

        'outer: for (i, &p) in Self::PRIMES.iter().enumerate() {
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

impl From<u8> for Denom {
    fn from(value: u8) -> Denom {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        if value <= 90 {
            return Denom::DECOMPOSED[value as usize - 1].clone();
        }
        Denom::decompose_u8(value)
    }
}

impl From<u16> for Denom {
    fn from(value: u16) -> Denom {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        if value <= 90 {
            return Denom::DECOMPOSED[value as usize - 1].clone();
        }
        Denom::decompose_u16(value)
    }
}

impl From<u32> for Denom {
    fn from(value: u32) -> Denom {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        if value <= 90 {
            return Denom::DECOMPOSED[value as usize - 1].clone();
        }
        Denom::decompose_u32(value)
    }
}

impl From<u64> for Denom {
    fn from(value: u64) -> Denom {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        if value <= 90 {
            return Denom::DECOMPOSED[value as usize - 1].clone();
        }
        Denom::decompose_u64(value)
    }
}

impl From<u128> for Denom {
    fn from(value: u128) -> Denom {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        if value <= 90 {
            return Denom::DECOMPOSED[value as usize - 1].clone();
        }
        Denom::decompose_u128(value)
    }
}

impl From<usize> for Denom {
    fn from(value: usize) -> Denom {
        if value == 0 {
            panic!("Attempted to create a denominator of zero");
        }
        if value <= 90 {
            return Denom::DECOMPOSED[value - 1].clone();
        }
        Denom::decompose_usize(value)
    }
}

impl From<BigUint> for Denom {
    fn from(value: BigUint) -> Denom {
        if value.is_zero() {
            panic!("Attempted to create a denominator of zero");
        }
        if value <= BigUint::from(90_usize) {
            return Denom::DECOMPOSED[TryInto::<usize>::try_into(value).unwrap() - 1].clone();
        }
        Denom::decompose_now(value)
    }
}

impl From<&BigUint> for Denom {
    fn from(value: &BigUint) -> Denom {
        if value.is_zero() {
            panic!("Attempted to create a denominator of zero");
        }
        if *value <= BigUint::from(90_usize) {
            return Denom::DECOMPOSED[TryInto::<usize>::try_into(value).unwrap() - 1].clone();
        }
        Denom::decompose_now(value.clone())
    }
}

impl From<Denom> for BigUint {
    fn from(value: Denom) -> BigUint {
        let mut result = match value.remainder {
            Some(x) => x,
            None => BigUint::ONE,
        };
        let mut tmp = 1_usize;
        for (i, &count) in value.primes.iter().enumerate() {
            let p = Denom::PRIMES[i];
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

impl From<&Denom> for BigUint {
    fn from(value: &Denom) -> BigUint {
        let mut result = match &value.remainder {
            Some(x) => x.clone(),
            None => BigUint::ONE,
        };
        let mut tmp = 1_usize;
        for (i, &count) in value.primes.iter().enumerate() {
            let p = Denom::PRIMES[i];
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

impl One for Denom {
    fn one() -> Self {
        Self::ONE
    }
}

impl num_traits::Pow<u32> for Denom {
    type Output = Self;

    fn pow(mut self, rhs: u32) -> Self {
        if rhs == 0 {
            return Denom::ONE;
        }
        let primes = Denom::pow_primes(&self.primes, rhs, &mut self.remainder);
        Denom {
            primes,
            remainder: self.remainder,
        }
    }
}

impl num_traits::Pow<u32> for &Denom {
    type Output = Denom;

    fn pow(self, rhs: u32) -> Denom {
        if rhs == 0 {
            return Denom::ONE;
        }
        let mut remainder = self.remainder.clone();
        let primes = Denom::pow_primes(&self.primes, rhs, &mut remainder);
        Denom { primes, remainder }
    }
}

impl Mul for Denom {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut remainder = match (self.remainder, rhs.remainder) {
            (None, None) => None,
            (None, Some(r)) => Some(r),
            (Some(l), None) => Some(l),
            (Some(l), Some(r)) => Some(l * r),
        };
        let primes = Denom::mul_primes(&self.primes, &rhs.primes, &mut remainder);
        Denom { primes, remainder }
    }
}

impl Mul<&Denom> for Denom {
    type Output = Self;

    fn mul(self, rhs: &Denom) -> Self {
        let mut remainder = match (self.remainder, &rhs.remainder) {
            (None, None) => None,
            (None, Some(r)) => Some(r.clone()),
            (Some(l), None) => Some(l),
            (Some(l), Some(r)) => Some(l * r),
        };
        let primes = Denom::mul_primes(&self.primes, &rhs.primes, &mut remainder);
        Denom { primes, remainder }
    }
}

impl Mul for &Denom {
    type Output = Denom;

    fn mul(self, rhs: Self) -> Denom {
        let mut remainder = match (&self.remainder, &rhs.remainder) {
            (None, None) => None,
            (None, Some(r)) => Some(r.clone()),
            (Some(l), None) => Some(l.clone()),
            (Some(l), Some(r)) => Some(l * r),
        };
        let primes = Denom::mul_primes(&self.primes, &rhs.primes, &mut remainder);
        Denom { primes, remainder }
    }
}

impl Mul<Denom> for &Denom {
    type Output = Denom;

    fn mul(self, rhs: Denom) -> Denom {
        let mut remainder = match (&self.remainder, rhs.remainder) {
            (None, None) => None,
            (None, Some(r)) => Some(r),
            (Some(l), None) => Some(l.clone()),
            (Some(l), Some(r)) => Some(l * r),
        };
        let primes = Denom::mul_primes(&self.primes, &rhs.primes, &mut remainder);
        Denom { primes, remainder }
    }
}

impl MulAssign for Denom {
    fn mul_assign(&mut self, rhs: Self) {
        match (&mut self.remainder, rhs.remainder) {
            (_, None) => (),
            (None, Some(r)) => self.remainder = Some(r),
            (Some(l), Some(r)) => *l *= r,
        };
        self.primes = Denom::mul_primes(&self.primes, &rhs.primes, &mut self.remainder);
    }
}

impl MulAssign<&Denom> for Denom {
    fn mul_assign(&mut self, rhs: &Denom) {
        match (&mut self.remainder, &rhs.remainder) {
            (_, None) => (),
            (None, Some(r)) => self.remainder = Some(r.clone()),
            (Some(l), Some(r)) => *l *= r,
        };
        self.primes = Denom::mul_primes(&self.primes, &rhs.primes, &mut self.remainder);
    }
}

impl Div for Denom {
    type Output = Self;

    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div(self, rhs: Self) -> Self {
        let (
            Denom {
                mut primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = Denom::div_primes(&self.primes, &rhs.primes, rhs.remainder);

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

        Denom::decompose_mask(&mut remainder, &mut primes, mask);
        Denom { primes, remainder }
    }
}

impl Div<&Denom> for Denom {
    type Output = Self;

    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div(self, rhs: &Denom) -> Self {
        let (
            Denom {
                mut primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = Denom::div_primes(&self.primes, &rhs.primes, rhs.remainder.clone());

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

        Denom::decompose_mask(&mut remainder, &mut primes, mask);
        Denom { primes, remainder }
    }
}

impl Div for &Denom {
    type Output = Denom;

    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div(self, rhs: Self) -> Denom {
        let (
            Denom {
                mut primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = Denom::div_primes(&self.primes, &rhs.primes, rhs.remainder.clone());

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

        Denom::decompose_mask(&mut remainder, &mut primes, mask);
        Denom { primes, remainder }
    }
}

impl Div<Denom> for &Denom {
    type Output = Denom;

    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div(self, rhs: Denom) -> Denom {
        let (
            Denom {
                mut primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = Denom::div_primes(&self.primes, &rhs.primes, rhs.remainder);

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

        Denom::decompose_mask(&mut remainder, &mut primes, mask);
        Denom { primes, remainder }
    }
}

impl DivAssign for Denom {
    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div_assign(&mut self, rhs: Self) {
        let (
            Denom {
                primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = Denom::div_primes(&self.primes, &rhs.primes, rhs.remainder);
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

        Denom::decompose_mask(&mut self.remainder, &mut self.primes, mask);
    }
}

impl DivAssign<&Denom> for Denom {
    /// Divides this denominator by the other one.
    ///
    /// This function panics if the other one doesn't divide this one.
    fn div_assign(&mut self, rhs: &Denom) {
        let (
            Denom {
                primes,
                remainder: rhs_remainder,
            },
            mask,
        ) = Denom::div_primes(&self.primes, &rhs.primes, rhs.remainder.clone());
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

        Denom::decompose_mask(&mut self.remainder, &mut self.primes, mask);
    }
}

impl Mul<Denom> for BigInt {
    type Output = Self;

    fn mul(mut self, rhs: Denom) -> Self {
        let mut tmp = 1_usize;
        for (i, &p) in Denom::PRIMES.iter().enumerate() {
            let count = rhs.primes[i];
            if count != 0 {
                Denom::accum_pow(&mut self, &mut tmp, p, count);
            }
        }

        self *= tmp;
        if let Some(remainder) = rhs.remainder {
            self *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        self
    }
}

impl Mul<&Denom> for BigInt {
    type Output = Self;

    fn mul(mut self, rhs: &Denom) -> Self {
        let mut tmp = 1_usize;
        for (i, &p) in Denom::PRIMES.iter().enumerate() {
            let count = rhs.primes[i];
            if count != 0 {
                Denom::accum_pow(&mut self, &mut tmp, p, count);
            }
        }

        self *= tmp;
        if let Some(remainder) = &rhs.remainder {
            self *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        self
    }
}

impl Mul<Denom> for &BigInt {
    type Output = BigInt;

    fn mul(self, rhs: Denom) -> BigInt {
        let mut this = self.clone();
        let mut tmp = 1_usize;
        for (i, &p) in Denom::PRIMES.iter().enumerate() {
            let count = rhs.primes[i];
            if count != 0 {
                Denom::accum_pow(&mut this, &mut tmp, p, count);
            }
        }

        this *= tmp;
        if let Some(remainder) = rhs.remainder {
            this *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        this
    }
}

impl Mul<&Denom> for &BigInt {
    type Output = BigInt;

    fn mul(self, rhs: &Denom) -> BigInt {
        let mut this = self.clone();
        let mut tmp = 1_usize;
        for (i, &p) in Denom::PRIMES.iter().enumerate() {
            let count = rhs.primes[i];
            if count != 0 {
                Denom::accum_pow(&mut this, &mut tmp, p, count);
            }
        }

        this *= tmp;
        if let Some(remainder) = &rhs.remainder {
            this *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        this
    }
}

impl Mul<BigInt> for Denom {
    type Output = BigInt;

    fn mul(self, mut rhs: BigInt) -> BigInt {
        let mut tmp = 1_usize;
        for (i, &p) in Denom::PRIMES.iter().enumerate() {
            let count = self.primes[i];
            if count != 0 {
                Denom::accum_pow(&mut rhs, &mut tmp, p, count);
            }
        }

        rhs *= tmp;
        if let Some(remainder) = self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        rhs
    }
}

impl Mul<&BigInt> for Denom {
    type Output = BigInt;

    fn mul(self, rhs: &BigInt) -> BigInt {
        let mut rhs = rhs.clone();
        let mut tmp = 1_usize;
        for (i, &p) in Denom::PRIMES.iter().enumerate() {
            let count = self.primes[i];
            if count != 0 {
                Denom::accum_pow(&mut rhs, &mut tmp, p, count);
            }
        }

        rhs *= tmp;
        if let Some(remainder) = self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder);
        }
        rhs
    }
}

impl Mul<BigInt> for &Denom {
    type Output = BigInt;

    fn mul(self, mut rhs: BigInt) -> BigInt {
        let mut tmp = 1_usize;
        for (i, &p) in Denom::PRIMES.iter().enumerate() {
            let count = self.primes[i];
            if count != 0 {
                Denom::accum_pow(&mut rhs, &mut tmp, p, count);
            }
        }

        rhs *= tmp;
        if let Some(remainder) = &self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        rhs
    }
}

impl Mul<&BigInt> for &Denom {
    type Output = BigInt;

    fn mul(self, rhs: &BigInt) -> BigInt {
        let mut rhs = rhs.clone();
        let mut tmp = 1_usize;
        for (i, &p) in Denom::PRIMES.iter().enumerate() {
            let count = self.primes[i];
            if count != 0 {
                Denom::accum_pow(&mut rhs, &mut tmp, p, count);
            }
        }

        rhs *= tmp;
        if let Some(remainder) = &self.remainder {
            rhs *= BigInt::from_biguint(Sign::Plus, remainder.clone());
        }
        rhs
    }
}

impl MulAssign<Denom> for BigInt {
    fn mul_assign(&mut self, rhs: Denom) {
        let mut tmp = 1_usize;
        for (i, &p) in Denom::PRIMES.iter().enumerate() {
            let count = rhs.primes[i];
            if count != 0 {
                Denom::accum_pow(self, &mut tmp, p, count);
            }
        }

        *self *= tmp;
        if let Some(remainder) = rhs.remainder {
            *self *= BigInt::from_biguint(Sign::Plus, remainder);
        }
    }
}

impl MulAssign<&Denom> for BigInt {
    fn mul_assign(&mut self, rhs: &Denom) {
        let mut tmp = 1_usize;
        for (i, &p) in Denom::PRIMES.iter().enumerate() {
            let count = rhs.primes[i];
            if count != 0 {
                Denom::accum_pow(self, &mut tmp, p, count);
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
    fn test_decompose_small() {
        for (i, x) in Denom::DECOMPOSED.iter().enumerate() {
            assert_eq!(x, &Denom::decompose_small(i + 1));
            assert_eq!(x, &Denom::from(BigUint::from(i + 1)));
            assert_eq!(x, &Denom::decompose_now(BigUint::from(i + 1)));
        }
    }

    #[test]
    #[should_panic(expected = "Failed to decompose small integer into small prime factors.")]
    fn test_decompose_small_out_of_range() {
        assert_eq!(
            Denom::decompose_small(97).to_biguint(),
            BigUint::from(97_usize)
        );
    }

    #[test]
    fn test_decompose_is_correct() {
        for i in 1_usize..=1000 {
            let bigi = BigUint::from(i);
            let x = Denom::from(&bigi);
            let mut recomposed = x.remainder.unwrap_or_else(BigUint::one);
            for (i, &prime) in Denom::PRIMES.iter().enumerate() {
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
            Denom::from(BigUint::from(128_usize)),
            Denom {
                primes: [
                    7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                ],
                remainder: None,
            }
        );
        assert_eq!(
            Denom::from(BigUint::from(89_usize)),
            Denom {
                primes: [
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
                ],
                remainder: None,
            }
        );
        assert_eq!(
            Denom::from(BigUint::from(97_usize)),
            Denom {
                primes: [
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                ],
                remainder: Some(BigUint::from(97_usize)),
            }
        );
        assert_eq!(
            Denom::from(BigUint::from(97000_usize)),
            Denom {
                primes: [
                    3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                ],
                remainder: Some(BigUint::from(97_usize)),
            }
        );
    }

    #[test]
    fn test_decompose_prime_powers() {
        for (i, &p) in Denom::PRIMES.iter().enumerate() {
            let p = BigUint::from(p);
            for power in 1..=255 {
                assert_eq!(
                    Denom::from(p.pow(power as u32)),
                    Denom {
                        primes: std::array::from_fn(|j| if i == j { power } else { 0 }),
                        remainder: None,
                    }
                );
            }
            for power in 1..=64 {
                assert_eq!(
                    Denom::from(p.pow(255 + power)),
                    Denom {
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
            let denom_a = Denom::from(p.pow(a));
            for b in 1..=256 {
                let denom_b = Denom::from(p.pow(b));
                let denom_ab = Denom::from(p.pow(a + b));
                assert_eq!(&denom_a * denom_b, denom_ab);
            }
        }
    }

    #[test]
    fn test_div_prime_powers() {
        let p = BigUint::from(2u32);
        for a in 1..=256 {
            let denom_a = Denom::from(p.pow(a));
            for b in 1..=256 {
                let denom_b = Denom::from(p.pow(b));
                let denom_ab = Denom::from(p.pow(a + b));
                assert_eq!(denom_ab / denom_b, denom_a);
            }
        }
    }

    #[test]
    fn test_to_biguint() {
        for i in 1_usize..=1000 {
            let bigi = BigUint::from(i);
            let x = Denom::from(&bigi);
            assert_eq!(x.to_biguint(), bigi);
        }
    }

    #[test]
    fn test_product() {
        let values = (100..200)
            .map(|i: usize| Denom::from(BigUint::from(i)))
            .collect::<Vec<_>>();
        for (i, x) in values.iter().enumerate().map(|(i, x)| (i + 100, x)) {
            for (j, y) in values.iter().enumerate().map(|(j, y)| (j + 100, y)) {
                let z = x * y;
                assert_eq!(z, Denom::from(BigUint::from(i * j)));
                for k in 0..Denom::NUM_PRIMES {
                    assert_eq!(z.primes[k], x.primes[k] + y.primes[k]);
                }
            }
        }
    }

    #[test]
    fn test_normalize() {
        let values = (100..200)
            .map(|i: usize| Denom::from(BigUint::from(i)))
            .collect::<Vec<_>>();
        for x in &values {
            for y in &values {
                let mut xnum = BigInt::one();
                let mut ynum = BigInt::one();
                let lcm = Denom::normalize(&mut xnum, &mut ynum, x, y);
                let lcm_bigint = lcm.to_biguint();
                let xnum = xnum.to_biguint().unwrap();
                let ynum = ynum.to_biguint().unwrap();

                assert_eq!(xnum * x.to_biguint(), lcm_bigint);
                assert_eq!(ynum * y.to_biguint(), lcm_bigint);
                for k in 0..Denom::NUM_PRIMES {
                    assert_eq!(lcm.primes[k], std::cmp::max(x.primes[k], y.primes[k]));
                }
            }
        }
    }

    #[test]
    fn test_gcd_reduce() {
        let mut num = BigInt::from(-3 * 97);
        let mut denom = Denom::from(BigUint::from(3u32 * 5 * 97));
        denom.gcd_reduce(&mut num);
        assert_eq!(num, BigInt::from(-1));
        assert_eq!(denom, Denom::from(BigUint::from(5u32)));
    }
}
