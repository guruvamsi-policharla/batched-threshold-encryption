use crate::utils::transpose;
use ark_ec::{pairing::Pairing, scalar_mul::fixed_base::FixedBase, Group};
use ark_ff::PrimeField;
use ark_poly::{domain::EvaluationDomain, Radix2EvaluationDomain};
use ark_std::{end_timer, rand::RngCore, start_timer, One, UniformRand, Zero};
use std::{iter, marker::PhantomData, vec};

pub struct CRS<E: Pairing> {
    pub powers_of_g: Vec<E::G1>,
    pub h: E::G2,
}

pub struct Dealer<E: Pairing> {
    batch_size: usize,
    n: usize,
    long_term_secret: E::ScalarField, // L_0(1) 
    _engine: PhantomData<E>,
}

impl<E> Dealer<E>
where
    E: Pairing,
{
    pub fn new(batch_size: usize, n: usize) -> Self {
        Self {
            batch_size,
            n,
            long_term_secret: E::ScalarField::zero(),
            _engine: PhantomData,
        }
    }

    pub fn setup<R: RngCore>(
        &mut self,
        rng: &mut R,
    ) -> (CRS<E>, Vec<Vec<E::ScalarField>>) {
        // Sample tau and compute its powers
        let tau = E::ScalarField::rand(rng);
        let powers_of_tau: Vec<E::ScalarField> =
            iter::successors(Some(E::ScalarField::one()), |p| Some(*p * tau))
                .take(self.batch_size + 1)
                .collect();

        // Generators
        let g = E::G1::generator();
        let h = E::G2::generator();

        // Compute powers of g
        let window_size = FixedBase::get_mul_window_size(self.batch_size + 1);
        let scalar_bits = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let g_time = start_timer!(|| "Generating powers of G");
        let g_table = FixedBase::get_window_table(scalar_bits, window_size, g);
        let powers_of_g =
            FixedBase::msm::<E::G1>(scalar_bits, window_size, &g_table, &powers_of_tau);
        end_timer!(g_time);
        
        // Secret share the lagrange coefficients for the domain {tau, omega, omega^2, ... omega^{batch_size}} at the evaluation point 1
        let tx_domain = Radix2EvaluationDomain::<E::ScalarField>::new(self.batch_size+1).unwrap();
        debug_assert!(tx_domain.size() == self.batch_size+1, "Could not find somain of batch_size+1 {}!", self.batch_size+1);
        
        let omega = tx_domain.group_gen;
        let mut lag_domain: Vec<E::ScalarField> = vec![E::ScalarField::one(); self.batch_size+1];

        lag_domain[0] = E::ScalarField::one();
        for i in 1..self.batch_size+1 {
            lag_domain[i] = lag_domain[i-1] * omega;
        }
        lag_domain[0] = tau; //lag_domain now contains {tau, omega, omega^2, ... omega^{batch_size}}

        // compute lagrange coefficients at these points
        // takes O(n^2) time but only done once
        // can be improved to O(nlogn) using FFT on a good domain first and then swapping out one point in each coefficient
        // and one entire coefficient
        let mut lag_coeffs: Vec<E::ScalarField> = vec![E::ScalarField::one(); self.batch_size+1];
        for i in 0..self.batch_size+1 {
            let mut num = E::ScalarField::one();
            let mut den = E::ScalarField::one();
            for j in 0..self.batch_size+1 {
                if i != j {
                    num *= E::ScalarField::one() - lag_domain[j];
                    den *= lag_domain[i] - lag_domain[j];
                }
            }
            lag_coeffs[i] = num / den;
        }

        // save the long_term_secret L_0(1)
        self.long_term_secret = lag_coeffs[0];

        // secret share the lagrange coefficients
        let mut lag_shares: Vec<Vec<E::ScalarField>> = Vec::new();
        let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(self.n).unwrap();
        for i in 1..=self.batch_size {
            // secret share i-th coefficient
            let mut coeffs = vec![E::ScalarField::zero(); self.n];
            coeffs[0] = powers_of_tau[i];
            for j in 1..self.n {
                coeffs[j] = E::ScalarField::rand(rng);
            }

            // use fft to compute shares
            let evals = share_domain.fft(&coeffs);
            lag_shares.push(evals);
        }

        let crs = CRS::<E> { powers_of_g, h };

        (crs, transpose(lag_shares))
    }

    pub fn epoch_setup<R: RngCore>(
        &self,
        rng: &mut R,
    ) -> (E::G1, Vec<E::ScalarField>){
        // secret share alpha*long_term_secret for a random alpha
        // publish com = g^alpha

        let alpha = E::ScalarField::rand(rng);
        let com = E::G1::generator() * alpha;

        let mut coeffs = vec![E::ScalarField::zero(); self.n];
        coeffs[0] = alpha * self.long_term_secret;
        for j in 1..self.n {
            coeffs[j] = E::ScalarField::rand(rng);
        }

        // use fft to compute shares
        let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(self.n).unwrap();
        let evals = share_domain.fft(&coeffs);

        (com, evals)
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use ark_bls12_381::Bls12_381;

    #[test]
    fn test_dealer() {
        let mut rng = ark_std::test_rng();
        let batch_size = (1 << 5) - 1;
        let n = 1 << 4;

        let mut dealer = Dealer::<Bls12_381>::new(batch_size, n);
        let (crs, lag_shares) = dealer.setup(&mut rng);
        let dealer = dealer;

        let (_com, evals) = dealer.epoch_setup(&mut rng);

        assert_eq!(lag_shares.len(), n);
        assert_eq!(lag_shares[0].len(), batch_size);
        assert_eq!(crs.powers_of_g.len(), batch_size + 1);
        assert_eq!(evals.len(), n);
    }
}