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

include!(concat!(env!("OUT_DIR"), "/primes.rs"));

pub fn known_odd_prime_factor_indices(x: usize) -> Option<impl Iterator<Item = (u16, u16)>> {
    ODD_FACTOR_INDICES
        .get(x / 2)
        .map(|array| array.iter().copied().filter(|(p, _)| *p != 0))
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::util::OddDivider;

    #[test]
    fn test_dividers_u16() {
        for (i, &p) in ODD_PRIMES.iter().enumerate() {
            let divider = OddDivider {
                divisor: p,
                multiplier: ODD_PRIME_DIVIDERS_U16[i],
                shift: p.ilog2(),
            };
            for x in 0..=u16::MAX {
                let (quo, rem) = divider.div_rem(x);
                assert_eq!(x / p, quo);
                assert_eq!(x % p, rem);
            }
        }
    }

    #[test]
    fn test_dividers_u32() {
        for (i, &p) in ODD_PRIMES.iter().enumerate() {
            let p = p as u32;
            let divider = OddDivider {
                divisor: p,
                multiplier: ODD_PRIME_DIVIDERS_U32[i],
                shift: p.ilog2(),
            };
            for x in (0..=u16::MAX as u32)
                .chain((0..=u16::MAX as u32).map(|x| u32::MAX - x))
                .chain((0..32).map(|i| 1 << i))
                .chain((0..32).map(|i| !(1 << i)))
                .chain((0..32).map(|i| (1 << i) - 1))
            {
                let (quo, rem) = divider.div_rem(x);
                assert_eq!(x / p, quo);
                assert_eq!(x % p, rem);
            }
        }
    }

    #[test]
    fn test_dividers_u64() {
        for (i, &p) in ODD_PRIMES.iter().enumerate() {
            let p = p as u64;
            let divider = OddDivider {
                divisor: p,
                multiplier: ODD_PRIME_DIVIDERS_U64[i],
                shift: p.ilog2(),
            };
            for x in (0..=u16::MAX as u64)
                .chain((0..=u16::MAX as u64).map(|x| u64::MAX - x))
                .chain((0..64).map(|i| 1 << i))
                .chain((0..64).map(|i| !(1 << i)))
                .chain((0..64).map(|i| (1 << i) - 1))
            {
                let (quo, rem) = divider.div_rem(x);
                assert_eq!(x / p, quo);
                assert_eq!(x % p, rem);
            }
        }
    }
}
