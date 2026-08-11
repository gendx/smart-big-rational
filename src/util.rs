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

use std::ops::{Div, Rem, Shr};

#[derive(Debug, Default, Clone, Copy)]
pub struct OddDivider<T> {
    /// Original divisor.
    pub divisor: T,
    /// Magic multiplier to compute divisions by this divisor.
    pub multiplier: T,
    /// Shift to compute divisions by this divisor.
    pub shift: u32,
}

impl<T> OddDivider<T>
where
    T: Copy + Shr<u32, Output = T> + Div<Output = T> + Rem<Output = T> + Arithmetic,
{
    // Not suitable for powers of 2 (as the shift is adjusted differently).
    #[inline(always)]
    fn div_non_power_of_two(&self, x: T) -> T {
        // See https://rubenvannieuwpoort.nl/posts/division-by-constant-unsigned-integers:
        //  multiplier = m - 2^N
        //  hi = (x*m - x*2^N) >> N = x*m / 2^N - x
        //  y = (2x - x*m / 2^N) / 2 + x*m / 2^N - x
        //    = x*m / 2^(N+1)
        //  y >> shift = x*m / 2^(N+shift+1)
        let (_, hi) = x.widening_mul(self.multiplier);
        let y = ((x.wrapping_sub(hi)) >> 1).wrapping_add(hi);
        y >> self.shift
    }

    #[inline(always)]
    pub fn div_rem(&self, x: T) -> (T, T) {
        let q = self.div_non_power_of_two(x);
        let r = x.wrapping_sub(q.wrapping_mul(self.divisor));
        (q, r)
    }
}

pub trait Arithmetic: Sized {
    fn wrapping_add(self, other: Self) -> Self;

    fn wrapping_sub(self, other: Self) -> Self;

    fn wrapping_mul(self, other: Self) -> Self;

    /// Returns (low, high) of multiplying `self` by `other`.
    fn widening_mul(self, other: Self) -> (Self, Self);
}

impl Arithmetic for u16 {
    #[inline(always)]
    fn wrapping_add(self, other: Self) -> Self {
        self.wrapping_add(other)
    }

    #[inline(always)]
    fn wrapping_sub(self, other: Self) -> Self {
        self.wrapping_sub(other)
    }

    #[inline(always)]
    fn wrapping_mul(self, other: Self) -> Self {
        self.wrapping_mul(other)
    }

    #[inline(always)]
    fn widening_mul(self, other: Self) -> (Self, Self) {
        self.carrying_mul(other, 0)
    }
}

impl Arithmetic for u32 {
    #[inline(always)]
    fn wrapping_add(self, other: Self) -> Self {
        self.wrapping_add(other)
    }

    #[inline(always)]
    fn wrapping_sub(self, other: Self) -> Self {
        self.wrapping_sub(other)
    }

    #[inline(always)]
    fn wrapping_mul(self, other: Self) -> Self {
        self.wrapping_mul(other)
    }

    #[inline(always)]
    fn widening_mul(self, other: Self) -> (Self, Self) {
        self.carrying_mul(other, 0)
    }
}

impl Arithmetic for u64 {
    #[inline(always)]
    fn wrapping_add(self, other: Self) -> Self {
        self.wrapping_add(other)
    }

    #[inline(always)]
    fn wrapping_sub(self, other: Self) -> Self {
        self.wrapping_sub(other)
    }

    #[inline(always)]
    fn wrapping_mul(self, other: Self) -> Self {
        self.wrapping_mul(other)
    }

    #[inline(always)]
    fn widening_mul(self, other: Self) -> (Self, Self) {
        self.carrying_mul(other, 0)
    }
}
