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
    sbox_degree: u64,
    constants: Vec<Vec<F>>,
    mat_internal: Vec<Vec<F>>,
}

impl<F, const WIDTH: usize> Poseidon2<F, WIDTH>
where
    F: Field,
{
    pub fn new(
        half_num_full_rounds: usize,
        num_partial_rounds: usize,
        sbox_degree: u64,
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
            self.full_sbox_layer(state);
            self.mul_external(state);
            *round_ctr += 1;
        }
    }

    fn partial_rounds(&self, state: &mut [F; WIDTH], round_ctr: &mut usize) {
        for _ in 0..self.num_partial_rounds {
            self.constant_layer(state, *round_ctr);
            self.partial_sbox_layer(state);
            self.mul_internal(state);
            *round_ctr += 1;
        }
    }

    fn full_sbox_layer(&self, state: &mut [F; WIDTH]) {
        for x in state.iter_mut() {
            *x = x.exp_u64(self.sbox_degree);
        }
    }

    fn partial_sbox_layer(&self, state: &mut [F; WIDTH]) {
        state[0] = state[0].exp_u64(self.sbox_degree);
    }

    fn mul_external(&self, state: &mut [F; WIDTH]) {
        match WIDTH {
            2 => {
                // Matrix circ(2, 1)
                let mut sum = state[0];
                sum.add_assign(state[1]);
                state[0].add_assign(sum);
                state[1].add_assign(sum);
            }
            3 => {
                // Matrix circ(2, 1, 1)
                let mut sum = state[0];
                sum.add_assign(state[1]);
                sum.add_assign(state[2]);
                state[0].add_assign(sum);
                state[1].add_assign(sum);
                state[2].add_assign(sum);
            }
            4 | 8 | 12 | 16 | 20 | 24 => {
                // Applying cheap 4x4 MDS matrix to each 4-element part of the state
                let t4 = WIDTH / 4;
                for i in 0..t4 {
                    let start_index = i * 4;
                    let mut t_0 = state[start_index];
                    t_0.add_assign(state[start_index + 1]);
                    let mut t_1 = state[start_index + 2];
                    t_1.add_assign(state[start_index + 3]);
                    let mut t_2 = state[start_index + 1];
                    t_2 = t_2 * 2; 
                    t_2.add_assign(t_1);
                    let mut t_3 = state[start_index + 3];
                    t_3.double_in_place();
                    t_3.add_assign(t_0);
                    let mut t_4 = t_1;
                    t_4.double_in_place();
                    t_4.double_in_place();
                    t_4.add_assign(t_3);
                    let mut t_5 = t_0;
                    t_5.double_in_place();
                    t_5.double_in_place();
                    t_5.add_assign(t_2);
                    let mut t_6 = t_3;
                    t_6.add_assign(t_5);
                    let mut t_7 = t_2;
                    t_7.add_assign(t_4);
                    state[start_index] = t_6;
                    state[start_index + 1] = t_5;
                    state[start_index + 2] = t_7;
                    state[start_index + 3] = t_4;
                }

                // Applying second cheap matrix
                let mut stored = [F::zero(); 4];
                for l in 0..4 {
                    stored[l] = state[l];
                    for j in 1..t4 {
                        stored[l].add_assign(&state[4 * j + l]);
                    }
                }
                for i in 0..state.len() {
                    state[i].add_assign(&stored[i % 4]);
                }
            }
            _ => {
                panic!()
            }
        }
    }

    fn mul_internal(&self, state: &mut [F; WIDTH]) {
        match WIDTH {
            2 => {
                // [2, 1]
                // [1, 3]
                let mut sum = state[0];
                sum.add_assign(state[1]);
                state[0].add_assign(sum);
                state[1].double_in_place();
                state[1].add_assign(sum);
            }
            3 => {
                // [2, 1, 1]
                // [1, 2, 1]
                // [1, 1, 3]
                let mut sum = state[0];
                sum.add_assign(state[1]);
                sum.add_assign(state[2]);
                state[0].add_assign(sum);
                state[1].add_assign(sum);
                state[2].double_in_place();
                state[2].add_assign(sum);
            }
            4 | 8 | 12 | 16 | 20 | 24 => {
                // Compute state sum
                let mut sum = state[0];
                state
                    .iter()
                    .skip(1)
                    .take(WIDTH-1)
                    .for_each(|el| sum.add_assign(*el));
                // Add sum + diag entry * element to each element
                for i in 0..state.len() {
                    state[i].mul_assign(self.mat_internal[i][i]);
                    state[i].add_assign(sum);
                }
            }
            _ => {
                panic!()
            }
        }
    }

    fn constant_layer(&self, state: &mut [F; WIDTH], round: usize) {
        state.iter().zip(self.constants[round].iter()).for_each(|(x, c)| {
            *x += *c;
        });
    }

    fn permute(&self, mut state: [F; WIDTH]) -> [F; WIDTH] {
        let mut round_ctr = 0;
        self.half_full_rounds(&mut state, &mut round_ctr);
        self.partial_rounds(&mut state, &mut round_ctr);
        self.half_full_rounds(&mut state, &mut round_ctr);
        state
    }
}

#[cfg(test)]
mod poseidon2_tests_goldilocks {
    use super::*;
    
    use crate::poseidon2::poseidon2_instance_goldilocks::{
        POSEIDON2_GOLDILOCKS_8_PARAMS,
        POSEIDON2_GOLDILOCKS_12_PARAMS,
        POSEIDON2_GOLDILOCKS_16_PARAMS,
        POSEIDON2_GOLDILOCKS_20_PARAMS,
    };
    use std::convert::TryFrom;

    type Scalar = FpGoldiLocks;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        use rand::{thread_rng, Rng}; 
        let instances = vec![
            Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS)
        ];
        for instance in instances {
            let t = instance.params.t;
            for _ in 0..TESTRUNS {
                let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

                let mut input2: Vec<Scalar>;
                loop {
                    input2 = (0..t).map(|_| F::rand(&mut rand::thread_rng())).collect();
                    if input1 != input2 {
                        break;
                    }
                }

                let perm1 = instance.permute(&input1);
                let perm2 = instance.permute(&input1);
                let perm3 = instance.permute(&input2);
                assert_eq!(perm1, perm2);
                assert_ne!(perm1, perm3);
            }
        }
    }
}
