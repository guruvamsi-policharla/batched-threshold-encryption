use crate::utils::transpose;
use ark_ec::{pairing::Pairing, scalar_mul::fixed_base::FixedBase, Group};
use ark_ff::{FftField, PrimeField};
use ark_poly::{domain::EvaluationDomain, Radix2EvaluationDomain};
use ark_std::{end_timer, rand::RngCore, start_timer, One, UniformRand, Zero};
use std::{iter, marker::PhantomData, vec};

pub struct CRS<E: Pairing> {
    pub powers_of_g: Vec<E::G1>,
    pub pk: E::G2,

    pub powers_of_top_tau: Vec<E::G1>,
}
pub struct Dealer<E: Pairing> {
    batch_size: usize,
    n: usize,
    pub long_term_secret: E::ScalarField, // L_0(1)
    pub tau: E::ScalarField,
    pub alpha: E::ScalarField,
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
            tau: E::ScalarField::zero(),
            alpha: E::ScalarField::zero(),
            _engine: PhantomData,
        }
    }

    pub fn setup<R: RngCore>(&mut self, rng: &mut R) -> (CRS<E>, Vec<Vec<E::ScalarField>>) {
        // Sample tau and compute its powers ==========================================================
        let tau = E::ScalarField::rand(rng);
        self.tau = tau;
        let powers_of_tau: Vec<E::ScalarField> =
            iter::successors(Some(E::ScalarField::one()), |p| Some(*p * tau))
                .take(self.batch_size + 1)
                .collect();

        // Generators
        let g = E::G1::generator();
        let h = E::G2::generator();

        // Compute powers of g
        let window_size = FixedBase::get_mul_window_size(self.batch_size + 1);
        let scalar_size = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let g_time = start_timer!(|| "Generating powers of G");
        let g_table = FixedBase::get_window_table(scalar_size, window_size, g);
        let powers_of_g =
            FixedBase::msm::<E::G1>(scalar_size, window_size, &g_table, &powers_of_tau);
        end_timer!(g_time);

        // Compute the Toeplitz matrix preprocessing ==================================================
        let mut top_tau = powers_of_tau.clone();
        top_tau.truncate(self.batch_size);
        top_tau.reverse();
        top_tau.resize(2 * self.batch_size, E::ScalarField::zero());

        let top_domain =
            Radix2EvaluationDomain::<E::ScalarField>::new(2 * self.batch_size).unwrap();
        let top_tau = top_domain.fft(&top_tau);

        // Compute powers of top_tau
        let window_size = FixedBase::get_mul_window_size(2 * self.batch_size);
        let scalar_size = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let top_tau_time = start_timer!(|| "Generating powers of top_tau");
        let top_tau_table = FixedBase::get_window_table(scalar_size, window_size, g);
        let powers_of_top_tau =
            FixedBase::msm::<E::G1>(scalar_size, window_size, &top_tau_table, &top_tau);
        end_timer!(top_tau_time);

        // Secret share the lagrange coefficients for the domain {1, omega, omega^2, ... omega^{batch_size-1}, tau} at the evaluation point gamma
        let tx_domain = Radix2EvaluationDomain::<E::ScalarField>::new(self.batch_size).unwrap();
        debug_assert!(
            tx_domain.size() == self.batch_size,
            "Could not find domain of batch_size {}!",
            self.batch_size
        );

        let mut lag_domain: Vec<E::ScalarField> = tx_domain.elements().collect();
        lag_domain.resize(self.batch_size + 1, E::ScalarField::zero());
        lag_domain[self.batch_size] = tau; //lag_domain now contains {1, omega, omega^2, ... omega^{batch_size-1}, tau}

        // compute lagrange coefficients at these points for evaluation at gamma
        // takes O(n^2) time but only done once
        // can be improved to O(nlogn) using FFT on a good domain first and then swapping out one point in each coefficient
        // and one entire coefficient
        let gamma = E::ScalarField::GENERATOR;
        let mut lag_coeffs: Vec<E::ScalarField> = vec![E::ScalarField::one(); self.batch_size + 1];
        for i in 0..self.batch_size + 1 {
            let mut num = E::ScalarField::one();
            let mut den = E::ScalarField::one();
            for j in 0..self.batch_size + 1 {
                if i != j {
                    num *= gamma - lag_domain[j];
                    den *= lag_domain[i] - lag_domain[j];
                }
            }
            lag_coeffs[i] = num / den;
        }

        // save the long_term_secret L_B(gamma)
        self.long_term_secret = lag_coeffs[self.batch_size];

        // secret share the lagrange coefficients
        let mut lag_shares: Vec<Vec<E::ScalarField>> = Vec::new();
        let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(self.n).unwrap();
        for i in 0..self.batch_size {
            // secret share i-th coefficient
            let mut coeffs = vec![E::ScalarField::zero(); self.n];
            coeffs[0] = powers_of_tau[i];
            // for j in 1..self.n {
            //     coeffs[j] = E::ScalarField::rand(rng);
            // }

            // use fft to compute shares
            let evals = share_domain.fft(&coeffs);
            lag_shares.push(evals);
        }

        let crs = CRS::<E> {
            powers_of_g,
            pk: h * tau,
            powers_of_top_tau,
        };

        (crs, transpose(lag_shares))
    }

    pub fn epoch_setup<R: RngCore>(&mut self, rng: &mut R) -> (E::G1, Vec<E::ScalarField>) {
        // secret share alpha*long_term_secret for a random alpha
        // publish com = g^alpha

        let alpha = E::ScalarField::rand(rng);
        self.alpha = alpha;
        let com = E::G1::generator() * alpha;

        let mut coeffs = vec![E::ScalarField::zero(); self.n];
        coeffs[0] = alpha * self.long_term_secret;
        // for j in 1..self.n {
        //     coeffs[j] = E::ScalarField::rand(rng);
        // }

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
        let batch_size = 1 << 5;
        let n = 1 << 4;

        let mut dealer = Dealer::<Bls12_381>::new(batch_size, n);
        let (crs, lag_shares) = dealer.setup(&mut rng);

        let (_com, evals) = dealer.epoch_setup(&mut rng);

        assert_eq!(lag_shares.len(), n);
        assert_eq!(lag_shares[0].len(), batch_size);
        assert_eq!(crs.powers_of_g.len(), batch_size + 1);
        assert_eq!(evals.len(), n);
    }
}
