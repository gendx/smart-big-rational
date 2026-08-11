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

use std::cmp::Ordering;
use std::env;
use std::fmt::Debug;
use std::fs::File;
use std::io::Write;
use std::iter::Peekable;
use std::ops::{Add, AddAssign, Shl, Shr};
use std::path::Path;

fn main() {
    let odd_primes = get_odd_primes();

    let out_dir = env::var_os("OUT_DIR").unwrap();
    let mut file = File::create(Path::new(&out_dir).join("primes.rs")).unwrap();

    writeln!(
        file,
        "/// All odd prime integers until 2^16. {} KiB.",
        odd_primes.len() >> 9,
    )
    .unwrap();
    writeln!(
        file,
        "pub static ODD_PRIMES: [u16; {}] = [",
        odd_primes.len(),
    )
    .unwrap();
    for line in odd_primes.chunks(10) {
        for (i, p) in line.iter().enumerate() {
            if i == 0 {
                write!(file, "    {p:#x},").unwrap();
            } else {
                write!(file, " {p:#x},").unwrap();
            }
        }
        writeln!(file).unwrap();
    }
    writeln!(file, "];").unwrap();

    writeln!(
        file,
        "/// Magic divider constants for odd prime integers until 2^16. {} KiB.",
        odd_primes.len() >> 9,
    )
    .unwrap();
    writeln!(
        file,
        "pub static ODD_PRIME_DIVIDERS_U16: [u16; {}] = [",
        odd_primes.len(),
    )
    .unwrap();
    for line in odd_primes.chunks(10) {
        for (i, &p) in line.iter().enumerate() {
            let q = odd_divider(p);
            if i == 0 {
                write!(file, "    {q:#x},").unwrap();
            } else {
                write!(file, " {q:#x},").unwrap();
            }
        }
        writeln!(file).unwrap();
    }
    writeln!(file, "];").unwrap();

    writeln!(
        file,
        "/// Magic divider constants for odd prime integers until 2^16. {} KiB.",
        odd_primes.len() >> 8,
    )
    .unwrap();
    writeln!(
        file,
        "pub static ODD_PRIME_DIVIDERS_U32: [u32; {}] = [",
        odd_primes.len(),
    )
    .unwrap();
    for line in odd_primes.chunks(10) {
        for (i, &p) in line.iter().enumerate() {
            let q = odd_divider(p as u32);
            if i == 0 {
                write!(file, "    {q:#x},").unwrap();
            } else {
                write!(file, " {q:#x},").unwrap();
            }
        }
        writeln!(file).unwrap();
    }
    writeln!(file, "];").unwrap();

    writeln!(
        file,
        "/// Magic divider constants for odd prime integers until 2^16. {} KiB.",
        odd_primes.len() >> 7,
    )
    .unwrap();
    writeln!(
        file,
        "pub static ODD_PRIME_DIVIDERS_U64: [u64; {}] = [",
        odd_primes.len(),
    )
    .unwrap();
    for line in odd_primes.chunks(5) {
        for (i, &p) in line.iter().enumerate() {
            let q = odd_divider(p as u64);
            if i == 0 {
                write!(file, "    {q:#x},").unwrap();
            } else {
                write!(file, " {q:#x},").unwrap();
            }
        }
        writeln!(file).unwrap();
    }
    writeln!(file, "];").unwrap();

    let mut prime_indices = [u16::MAX; 1 << 16];
    for (i, &p) in odd_primes.iter().enumerate() {
        prime_indices[p as usize] = i as u16 + 1;
    }

    let count_odd_factors = 4096;
    let factors = get_factors(count_odd_factors * 2);
    let factor_indices = get_factor_indices(&factors, &prime_indices);

    writeln!(
        file,
        "/// Prime indices (1 = 3, 2 = 5, 3 = 7, etc.) of factors of odd numbers until {}. {} KiB.",
        count_odd_factors * 2,
        count_odd_factors >> 6,
    )
    .unwrap();
    writeln!(
        file,
        "/// There are at most 4 distinct odd prime factors because 3*5*7*11*13 = 15015.",
    )
    .unwrap();
    writeln!(
        file,
        "static ODD_FACTOR_INDICES: [[(u16, u16); 4]; {count_odd_factors}] = [",
    )
    .unwrap();
    for line in factor_indices {
        write!(file, "    [").unwrap();
        assert_eq!(line.len(), 4);
        for (i, (p, count)) in line.iter().enumerate() {
            if i != 0 {
                write!(file, ", ").unwrap();
            }
            write!(file, "({p:#x}, {count})").unwrap();
        }
        writeln!(file, "],").unwrap();
    }
    writeln!(file, "];").unwrap();

    drop(file);
    println!("cargo::rerun-if-changed=build.rs");
}

fn get_odd_primes() -> Vec<u16> {
    let mut sieve = vec![true; 1 << 16];
    sieve[0] = false;
    sieve[1] = false;

    let mut odd_primes = Vec::new();
    for i in (3..u16::MAX).step_by(2) {
        if sieve[i as usize] {
            odd_primes.push(i);
            for j in (3..).step_by(2) {
                match i.checked_mul(j) {
                    None => break,
                    Some(product) => sieve[product as usize] = false,
                }
            }
        }
    }

    odd_primes
}

fn get_factor_indices(factors: &[[(u16, u16); 4]], prime_indices: &[u16]) -> Vec<[(u16, u16); 4]> {
    factors
        .iter()
        .map(|&array| {
            array.map(|(p, count)| {
                let index = if p == 0 { 0 } else { prime_indices[p as usize] };
                assert_ne!(index, u16::MAX, "index for {p}");
                (index, count)
            })
        })
        .collect()
}

fn get_factors(max: usize) -> Vec<[(u16, u16); 4]> {
    assert!(max <= (1 << 16));
    let mut table = vec![[(0, 0); 4]; max / 2];

    let last: u16 = (max - 1) as u16;
    for i in (3..=last).step_by(2) {
        if table[i as usize / 2][0].0 == 0 {
            table[i as usize / 2][0] = (i, 1);
        }
        for j in (3..=i).step_by(2) {
            match i.checked_mul(j) {
                Some(k) if k <= last => {
                    if table[k as usize / 2][0].0 == 0 {
                        table[k as usize / 2] = mul(&table[i as usize / 2], &table[j as usize / 2]);
                    }
                }
                _ => break,
            }
        }
    }

    table
}

fn mul(x: &[(u16, u16)], y: &[(u16, u16)]) -> [(u16, u16); 4] {
    let mut primes = Vec::new();
    for (p, (a, b)) in Zip(
        x.iter().copied().filter(|(p, _)| *p != 0).peekable(),
        y.iter().copied().filter(|(p, _)| *p != 0).peekable(),
    ) {
        primes.push((p, a + b));
    }
    assert!(primes.len() <= 4);
    for _ in primes.len()..4 {
        primes.push((0, 0));
    }
    primes.try_into().unwrap()
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

fn odd_divider<T>(divisor: T) -> T
where
    T: Copy
        + Eq
        + Ord
        + Debug
        + Add<Output = T>
        + AddAssign
        + Shl<u32, Output = T>
        + Shr<u32, Output = T>
        + Arithmetic,
{
    // Division by 0 is caught here.
    let log2 = divisor.checked_ilog2().unwrap();

    // See https://rubenvannieuwpoort.nl/posts/division-by-constant-unsigned-integers:
    //  multiplier = 2^(N+shift) / d
    //  rem = 2^(N+shift) % d
    let (multiplier, rem) = T::narrowing_div_rem((T::ZERO, T::ONE << log2), divisor);
    assert_ne!(multiplier, T::ZERO);
    assert!(rem > T::ZERO);
    assert!(rem < divisor);

    // At this point the highest bit of multiplier is set (because we divided
    // 2^(N+log2(d)) / d), therefore shifting discards it:
    //  multiplier = 2 * (2^(N+shift) / d) - 2^N
    assert_eq!(multiplier >> (T::BITS - 1), T::ONE);
    let mut multiplier = multiplier << 1;
    let twice_rem = rem << 1;
    // Use the remainder to adjust the multiplier to:
    //  multiplier = 2^(N+shift+1) / d - 2^N
    if twice_rem >= divisor || twice_rem < rem {
        multiplier += T::ONE;
    }

    // Lastly, we compute the ceiling of that:
    //  multiplier = 2^(N+shift+1) / d + 1 - 2^N
    // Because d isn't a power of 2 (and therefore doesn't divide 2^(N+shift+1)),
    // this gives:
    //  multiplier = ceil(2^(N+shift+1) / d) - 2^N
    multiplier + T::ONE
}

trait Arithmetic: Sized {
    const BITS: u32;
    const ZERO: Self;
    const ONE: Self;

    /// Returns the base-2 logarithm or `None` if `self` is zero.
    fn checked_ilog2(self) -> Option<u32>;

    /// Divides (num_lo, num_hi) by denom, returning (quotient, remainder),
    /// assuming that the quotient fits in Self.
    fn narrowing_div_rem(num: (Self, Self), denom: Self) -> (Self, Self);
}

impl Arithmetic for u16 {
    const BITS: u32 = u16::BITS;
    const ZERO: Self = 0;
    const ONE: Self = 1;

    #[inline(always)]
    fn checked_ilog2(self) -> Option<u32> {
        self.checked_ilog2()
    }

    #[inline(always)]
    fn narrowing_div_rem((num_lo, num_hi): (Self, Self), denom: Self) -> (Self, Self) {
        let a = ((num_hi as u32) << 16) | (num_lo as u32);
        let b = denom as u32;
        let quo = a / b;
        let rem = a.wrapping_sub(quo.wrapping_mul(b));
        (quo as u16, rem as u16)
    }
}

impl Arithmetic for u32 {
    const BITS: u32 = u32::BITS;
    const ZERO: Self = 0;
    const ONE: Self = 1;

    #[inline(always)]
    fn checked_ilog2(self) -> Option<u32> {
        self.checked_ilog2()
    }

    #[inline(always)]
    fn narrowing_div_rem((num_lo, num_hi): (Self, Self), denom: Self) -> (Self, Self) {
        let a = ((num_hi as u64) << 32) | (num_lo as u64);
        let b = denom as u64;
        let quo = a / b;
        let rem = a.wrapping_sub(quo.wrapping_mul(b));
        (quo as u32, rem as u32)
    }
}

impl Arithmetic for u64 {
    const BITS: u32 = u64::BITS;
    const ZERO: Self = 0;
    const ONE: Self = 1;

    #[inline(always)]
    fn checked_ilog2(self) -> Option<u32> {
        self.checked_ilog2()
    }

    #[inline(always)]
    fn narrowing_div_rem((num_lo, num_hi): (Self, Self), denom: Self) -> (Self, Self) {
        let a = ((num_hi as u128) << 64) | (num_lo as u128);
        let b = denom as u128;
        let quo = a / b;
        let rem = a.wrapping_sub(quo.wrapping_mul(b));
        (quo as u64, rem as u64)
    }
}
