// Copyright 2018 David Roundy. See the COPYRIGHT
// file at the top-level directory of this distribution and at
// https://rust-lang.org/COPYRIGHT.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! The [xoshiro128plus
//! RNG](http://xoshiro.di.unimi.it/xoroshiro128plus.c).

use std::num::Wrapping;
use std::{fmt};
use rand_core::{RngCore, SeedableRng, Error, impls, le};

/// The [xoshiro128plus
/// RNG](http://xoshiro.di.unimi.it/xoroshiro128plus.c).

#[derive(Clone)]
#[derive(Serialize,Deserialize)]
pub struct Xoroshiro128plusRng {
    s: [Wrapping<u64>; 2],
}

/// Our random number generator.
pub type Rng = Xoroshiro128plusRng;

// Custom Debug implementation that does not expose the internal state
impl fmt::Debug for Xoroshiro128plusRng {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Xoroshiro128plusRng {{}}")
    }
}

fn rotl(x: Wrapping<u64>, k: usize) -> Wrapping<u64> {
    (x << k) | (x >> (64-k))
}

impl RngCore for Xoroshiro128plusRng {
    #[inline]
    fn next_u32(&mut self) -> u32 {
        impls::next_u32_via_fill(self)
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        let s0 = self.s[0];
        let mut s1 = self.s[1];
        let result = s0 + s1;

        s1 ^= s0;
        self.s[0] = rotl(s0,24) ^ s1 ^ (s1 << 16); // a, b
        self.s[1] = rotl(s1, 37);
        result.0
    }

    #[inline]
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        impls::fill_bytes_via_next(self, dest)
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), Error> {
        Ok(self.fill_bytes(dest))
    }
}

impl Xoroshiro128plusRng {
    /// Seed this RNG using a u64.  This is not quite as trivial as
    /// you'd wish, because we need to ensure that the resulting state
    /// is not *all* zeros.
    pub fn from_u64(seed: u64) -> Self {
        let mut seed_u64 = [Wrapping(0u64); 2];
        // As recommended, we use splitmix64 to seed the generator if
        // we only have 64 bits of seed.
        let mut z = if seed == 0 {
            // Given a zero seed, we arbitrarily use all 1s, since
            // splitmix64 will return all zeros given a zero seed,
            // which is quite the opposite of what we want here.
            Wrapping(0xFFFFFFFFFFFFFFFF)
        } else {
            Wrapping(seed)
        };
	z = (z ^ (z >> 30)) * Wrapping(0xBF58476D1CE4E5B9);
	z = (z ^ (z >> 27)) * Wrapping(0x94D049BB133111EB);
        seed_u64[0] = z ^ (z >> 31);
	z = (z ^ (z >> 30)) * Wrapping(0xBF58476D1CE4E5B9);
	z = (z ^ (z >> 27)) * Wrapping(0x94D049BB133111EB);
        seed_u64[1] = z ^ (z >> 31);

        Xoroshiro128plusRng {
            s: seed_u64,
        }
    }
}

impl SeedableRng for Xoroshiro128plusRng {
    type Seed = [u8; 16];

    fn from_seed(seed: Self::Seed) -> Self {
        let mut seed_u64 = [0u64; 2];
        le::read_u64_into(&seed, &mut seed_u64);

        // As recommended, we use splitmix64 to seed the generator if
        // we only have 64 bits of seed.
        if seed_u64.iter().any(|&x| x == 0) {
            return Xoroshiro128plusRng::from_u64(seed_u64[0] | seed_u64[1]);
        }

        Xoroshiro128plusRng {
            s: [Wrapping(seed_u64[0]), Wrapping(seed_u64[1])],
        }
    }

    fn from_rng<R: RngCore>(mut rng: R) -> Result<Self, Error> {
        Ok(Xoroshiro128plusRng {
            s: [Wrapping(rng.next_u64()), Wrapping(rng.next_u64())],
        })
    }
}

#[cfg(test)]
mod tests {
    use rand_core::{RngCore, SeedableRng};
    use super::Xoroshiro128plusRng;

    #[test]
    fn test_xoshiro_zero_seed() {
        // Xoroshiro does not work with an all zero seed.
        // Assert it does not panic.
        let seed = [0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0];
        let mut rng = Xoroshiro128plusRng::from_seed(seed);
        let a = rng.next_u64();
        let b = rng.next_u64();
        assert!(a != 0);
        assert!(b != a);
    }

    #[test]
    fn test_xoshiro_clone() {
        let seed = [1,2,3,4, 5,5,7,8, 8,7,6,5, 4,3,2,1];
        let mut rng1 = Xoroshiro128plusRng::from_seed(seed);
        let mut rng2 = rng1.clone();
        for _ in 0..16 {
            assert_eq!(rng1.next_u64(), rng2.next_u64());
        }
    }
}
