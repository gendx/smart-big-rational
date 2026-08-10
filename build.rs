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

use std::env;
use std::fs::File;
use std::io::Write;
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
