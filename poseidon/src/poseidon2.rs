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
pub struct Poseidon2<F, const WIDTH: usize>
where
    F: Field,
{
    half_num_full_rounds: usize,
    num_partial_rounds: usize,
    sbox_degree: usize,
    constants: Vec<Vec<F>>,
    mat_internal: Vec<Vec<F>>,
}

impl<F, MDS, const WIDTH: usize> Poseidon2<F, WIDTH>
where
    F: Field,
{
    pub fn new(
        half_num_full_rounds: usize,
        num_partial_rounds: usize,
        sbox_degree: usize,
        constants: Vec<Vec<F>>,
        mat_internal: Vec<Vec<F>>,
    ) -> Self {
        let num_rounds = 2 * half_num_full_rounds + num_partial_rounds;
        assert_eq!(constants.len(), num_rounds);
        assert!(sbox_degree == 3 || sbox_degree == 5 || sbox_degree == 7 || sbox_degree == 11);
        Self {
            half_num_full_rounds,
            num_partial_rounds,
            sbox_degree,
            constants,
            mat_internal,
        }
    }

    fn half_full_rounds(&self, state: &mut [F; WIDTH], round_ctr: &mut usize) {
        for _ in 0..self.half_num_full_rounds {
            self.constant_layer(state, *round_ctr);
            Self::full_sbox_layer(state);
            self.mul_external(state);
            *round_ctr += 1;
        }
    }

    fn partial_rounds(&self, state: &mut [F; WIDTH], round_ctr: &mut usize) {
        for _ in 0..self.num_partial_rounds {
            self.constant_layer(state, *round_ctr);
            Self::partial_sbox_layer(state);
            self.mul_internal(state);
            *round_ctr += 1;
        }
    }

    fn full_sbox_layer(state: &mut [F; WIDTH]) {
        for x in state.iter_mut() {
            *x = x.exp_u64(ALPHA);
        }
    }

    fn partial_sbox_layer(state: &mut [F; WIDTH]) {
        state[0] = state[0].exp_u64(ALPHA);
    }

    fn mul_external(&self, state: &mut [F; WIDTH]) {
        let mut tmp = [F::ZERO; WIDTH];
        for i in 0..WIDTH {
            for j in 0..WIDTH {
                tmp[i] += self.mat_external[i][j] * state[j];
            }
        }
        *state = tmp;
    }

    fn mul_internal(&self, state: &mut [F; WIDTH]) {
        let mut tmp = [F::ZERO; WIDTH];
        for i in 0..WIDTH {
            for j in 0..WIDTH {
                tmp[i] += self.mat_internal[i][j] * state[j];
            }
        }
        *state = tmp;
    }

    fn constant_layer(&self, state: &mut [F; WIDTH], round: usize) {
        state.iter().zip(self.constants[round].iter()).for_each(|(x, c)| {
            *x += *c;
        });
    }
}

