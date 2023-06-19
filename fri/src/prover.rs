use crate::{FriConfig, FriProof};
use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Reverse;
use itertools::Itertools;
use p3_challenger::Challenger;
use p3_commit::{DirectMMCS, MMCS};
use p3_field::{AbstractField, ExtensionField, Field};
use p3_matrix::{Matrix, MatrixRows};
use p3_maybe_rayon::MaybeIntoParIter;
use p3_maybe_rayon::ParallelIterator;

pub(crate) fn prove<F, Challenge, M, MC, Chal>(
    codewords: &[M::ProverData],
    config: &FriConfig,
    challenger: &mut Chal,
) -> FriProof<F, Challenge, M, MC>
where
    F: Field,
    Challenge: ExtensionField<F>,
    M: MMCS<F>,
    MC: DirectMMCS<F>,
    Chal: Challenger<F>,
{
    let _commit_phase_commits = commit_phase::<F, Challenge, M, MC, Chal>(codewords, challenger);
    let query_indices: Vec<Challenge> = (0..config.num_queries)
        .map(|_| challenger.random_ext_element())
        .collect();
    let _query_proofs = query_indices.into_par_iter().map(|_query_index| todo!());
    todo!()
}

pub(crate) fn commit_phase<F, Challenge, M, MC, Chal>(
    codewords: &[M::ProverData],
    challenger: &mut Chal,
) -> Vec<MC::ProverData>
where
    F: Field,
    Challenge: ExtensionField<F>,
    M: MMCS<F>,
    MC: DirectMMCS<F>,
    Chal: Challenger<F>,
{
    let alpha: Challenge = challenger.random_ext_element();
    let matrices_by_desc_height = codewords
        .iter()
        .flat_map(|data| M::get_matrices(data))
        .sorted_by_key(|mat| Reverse(mat.height()))
        .group_by(|mat| mat.height());
    let mut matrices_by_desc_height = matrices_by_desc_height.into_iter();

    let (max_height, largest_matrices_iter) = matrices_by_desc_height.next().expect("No matrices?");
    let largest_matrices = largest_matrices_iter.collect_vec();
    let zero_vec = vec![<Challenge as AbstractField>::ZERO; max_height];
    let mut current = reduce_matrices(max_height, zero_vec, largest_matrices, alpha);
    let mut committed = vec![current.clone()];

    for (height, matrices) in matrices_by_desc_height {
        while current.len() < height {
            let beta = <Challenge as AbstractField>::ZERO; // TODO
            current = fold_even_odd(&current, beta);
        }
        committed.push(current.clone());
        current = reduce_matrices::<F, Challenge, M::Mat>(
            height,
            current.clone(),
            matrices.collect(),
            alpha,
        );
    }
    todo!()
}

/// Fold a polynomial `p(x) = p_even(x^2) + x p_odd(x^2)` into `p_even(x) + beta * p_odd(x)`.
fn fold_even_odd<F: Field>(_poly: &[F], _beta: F) -> Vec<F> {
    todo!()
}

#[allow(clippy::needless_pass_by_value)]
fn reduce_matrices<F, Challenge, Mat>(
    height: usize,
    init: Vec<Challenge>,
    matrices: Vec<&Mat>,
    alpha: Challenge,
) -> Vec<Challenge>
where
    F: Field,
    Challenge: ExtensionField<F>,
    Mat: for<'a> MatrixRows<'a, F>,
{
    (0..height)
        .map(|r| {
            let mut reduced = init[r];
            for mat in &matrices {
                for col in mat.row(r) {
                    reduced = reduced * alpha + *col;
                }
            }
            reduced
        })
        .collect()
}
