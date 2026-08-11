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

#![feature(test)]

extern crate test;

use ::test::Bencher;
use num_bigint::{BigInt, BigUint};
use smart_big_rational::Denom;
use std::hint::black_box;

macro_rules! benches {
    (
        $mod:ident,
        $denom:ty,
        $( $case:ident ,)*
    ) => {
        mod $mod {
            use super::*;
            use ::test::Bencher;

            $(
                #[bench]
                fn $case(b: &mut Bencher) {
                    $crate::$case::<$denom>(b);
                }
            )*
        }
    };
}

macro_rules! all_benches {
    (
        $mod:ident,
        $denom:ty
    ) => {
        benches!(
            $mod,
            $denom,
            decompose_large_prime_u016,
            decompose_large_prime_u032,
            decompose_large_prime_u064,
            decompose_large_prime_u128,
            decompose_large_prime_vbigint,
            decompose_small_factors_u008,
            decompose_small_factors_u016,
            decompose_small_factors_u032,
            decompose_small_factors_u064,
            decompose_small_factors_u128,
            decompose_small_factors_vbigint,
            decompose_large_factors_u016,
            decompose_large_factors_u032,
            decompose_large_factors_u064,
            decompose_large_factors_u128,
            decompose_large_factors_vbigint,
            decompose_all_u008,
            decompose_all_u016,
        );
    };
}

mod denom_array {
    use smart_big_rational::DenomArray;

    all_benches!(bench24, DenomArray<24>);
    all_benches!(bench6542, DenomArray<6542>);
}

mod denom_sparse {
    use smart_big_rational::DenomSparseU16;

    all_benches!(bench24, DenomSparseU16<24>);
    all_benches!(bench6542, DenomSparseU16<6542>);
}

fn decompose_large_prime_u016<D: Denom>(b: &mut Bencher) {
    let x: u16 = 0xfff1;
    b.iter(|| D::from(x));
}

fn decompose_large_prime_u032<D: Denom>(b: &mut Bencher) {
    let x: u32 = 0xffff_fffb;
    b.iter(|| D::from(x));
}

fn decompose_large_prime_u064<D: Denom>(b: &mut Bencher) {
    let x: u64 = 0xffff_ffff_ffff_ffc5;
    b.iter(|| D::from(x));
}

fn decompose_large_prime_u128<D: Denom>(b: &mut Bencher) {
    let x: u128 = 0xffff_ffff_ffff_ffff_ffff_ffff_ffff_ff61;
    b.iter(|| D::from(x));
}

fn decompose_large_prime_vbigint<D: Denom>(b: &mut Bencher) {
    let x: BigUint = BigUint::from_slice(&[
        0xffff_ff43,
        0xffff_ffff,
        0xffff_ffff,
        0xffff_ffff,
        0xffff_ffff,
        0xffff_ffff,
        0xffff_ffff,
        0xffff_ffff,
    ]);
    b.iter(|| D::from(&x));
}

fn decompose_small_factors_u008<D: Denom>(b: &mut Bencher) {
    let x: u8 = 2 * 3 * 5 * 7;
    b.iter(|| D::from(x));
}

fn decompose_small_factors_u016<D: Denom>(b: &mut Bencher) {
    let x: u16 = 2 * 3 * 5 * 7 * 11 * 13;
    b.iter(|| D::from(x));
}

fn decompose_small_factors_u032<D: Denom>(b: &mut Bencher) {
    let x: u32 = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23;
    b.iter(|| D::from(x));
}

fn decompose_small_factors_u064<D: Denom>(b: &mut Bencher) {
    let x: u64 = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 31 * 37 * 41 * 43 * 47 * 53;
    b.iter(|| D::from(x));
}

fn decompose_small_factors_u128<D: Denom>(b: &mut Bencher) {
    let x: u128 = 2
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
        * 101;
    b.iter(|| D::from(x));
}

fn decompose_small_factors_vbigint<D: Denom>(b: &mut Bencher) {
    let x: BigInt = BigInt::from(2)
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
        * 233;
    let x: BigUint = x.try_into().unwrap();
    b.iter(|| D::from(&x));
}

fn decompose_large_factors_u016<D: Denom>(b: &mut Bencher) {
    let x: u16 = 0xfb * 0xf1;
    b.iter(|| D::from(x));
}

fn decompose_large_factors_u032<D: Denom>(b: &mut Bencher) {
    let x: u32 = 0xfff1 * 0xffef;
    b.iter(|| D::from(x));
}

fn decompose_large_factors_u064<D: Denom>(b: &mut Bencher) {
    let x: u64 = 0xfff1 * 0xffef * 0xffd9 * 0xffc7;
    b.iter(|| D::from(x));
}

fn decompose_large_factors_u128<D: Denom>(b: &mut Bencher) {
    let x: u128 = 0xfff1 * 0xffef * 0xffd9 * 0xffc7 * 0xffa9 * 0xffa7 * 0xff9d * 0xff8f;
    b.iter(|| D::from(x));
}

fn decompose_large_factors_vbigint<D: Denom>(b: &mut Bencher) {
    let x = BigUint::from(0xfff1_u16)
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
        * 0xff49_u16;
    b.iter(|| D::from(&x));
}

fn decompose_all_u008<D: Denom>(b: &mut Bencher) {
    b.iter(|| {
        for x in 1..=u8::MAX {
            black_box(D::from(x));
        }
    });
}

fn decompose_all_u016<D: Denom>(b: &mut Bencher) {
    b.iter(|| {
        for x in 1..=u16::MAX {
            black_box(D::from(x));
        }
    });
}
