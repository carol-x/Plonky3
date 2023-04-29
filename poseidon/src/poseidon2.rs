use alloc::vec::Vec;
use p3_field::field::Field;
use p3_symmetric::permutation::{ArrayPermutation, CryptographicPermutation, MDSPermutation};
use p3_symmetric::sponge::PaddingFreeSponge;
use rand::distributions::Standard;
use rand::prelude::Distribution;
use rand::Rng;

// reference Poseidon2 paper: https://eprint.iacr.org/2023/323.pdf
// reference Poseidon2 code: https://github.com/HorizenLabs/poseidon2/blob/main/plain_implementations/src/poseidon2/poseidon2.rs

/// The Poseidon permutation 
