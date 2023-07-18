use crate::utils::transpose;
use ark_ec::{pairing::Pairing, scalar_mul::fixed_base::FixedBase, Group};
use ark_ff::PrimeField;
use ark_poly::{domain::EvaluationDomain, Radix2EvaluationDomain};
use ark_std::{end_timer, rand::RngCore, start_timer, One, UniformRand, Zero};
use std::{iter, marker::PhantomData};

pub struct CRS<E: Pairing> {
    pub powers_of_g: Vec<E::G1>,
    pub h: E::G2,
}

pub struct Dealer<E: Pairing> {
    _engine: PhantomData<E>,
}

impl<E> Dealer<E>
where
    E: Pairing,
{
    pub fn setup<R: RngCore>(
        batch_size: usize,
        n: usize,
        rng: &mut R,
    ) -> (CRS<E>, Vec<Vec<E::ScalarField>>) {
        // Sample tau and compute its powers
        let tau = E::ScalarField::rand(rng);
        let powers_of_tau: Vec<E::ScalarField> =
            iter::successors(Some(E::ScalarField::one()), |p| Some(*p * tau))
                .take(batch_size + 1)
                .collect();

        // Generators
        let g = E::G1::generator();
        let h = E::G2::generator();

        // Compute powers of g
        let window_size = FixedBase::get_mul_window_size(batch_size + 1);
        let scalar_bits = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let g_time = start_timer!(|| "Generating powers of G");
        let g_table = FixedBase::get_window_table(scalar_bits, window_size, g);
        let powers_of_g =
            FixedBase::msm::<E::G1>(scalar_bits, window_size, &g_table, &powers_of_tau);
        end_timer!(g_time);

        // Secrete share the powers of tau
        let mut tau_shares: Vec<Vec<E::ScalarField>> = Vec::new();
        let domain = Radix2EvaluationDomain::<E::ScalarField>::new(n).unwrap();

        for i in 1..=batch_size + 1 {
            // secret share i-th power of tau
            let mut coeffs = vec![E::ScalarField::zero(); n];
            coeffs[0] = powers_of_tau[i - 1];
            for j in 1..n {
                coeffs[j] = E::ScalarField::rand(rng);
            }

            // use fft to compute shares
            let evals = domain.fft(&coeffs);
            tau_shares.push(evals);
        }

        let crs = CRS::<E> { powers_of_g, h };

        (crs, transpose(tau_shares))
    }
}
