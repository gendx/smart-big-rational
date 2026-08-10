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
use std::fs::File;
use std::io::Write;
use std::iter::Peekable;
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
