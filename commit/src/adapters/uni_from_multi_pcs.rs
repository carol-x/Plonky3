use crate::pcs::MultivariatePCS;
use core::marker::PhantomData;
use p3_field::Field;

pub struct UniFromMultiPCS<F: Field, M: MultivariatePCS<F>> {
    _multi: M,
    _phantom_f: PhantomData<F>,
}

// impl<F: Field, M: MultivariatePCS<F>> UnivariatePCS<F> for UniFromMultiPCS<F> {}
