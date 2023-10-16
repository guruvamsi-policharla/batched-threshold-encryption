use crate::utils::{lagrange_coefficients, transpose};
use ark_ec::{pairing::Pairing, scalar_mul::fixed_base::FixedBase, Group};
use ark_ff::{FftField, PrimeField};
use ark_poly::{domain::EvaluationDomain, Radix2EvaluationDomain};
use ark_std::{rand::RngCore, One, UniformRand, Zero};
use std::{iter, marker::PhantomData, vec};
use ark_serialize::*;

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct CRS<E: Pairing> {
    pub powers_of_g: Vec<E::G1>,
    pub htau: E::G2,

    pub y: Vec<E::G1>,
}
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct Dealer<E: Pairing> {
    batch_size: usize,
    n: usize,
    pub long_term_secret: E::ScalarField, // L_0(1)
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

    pub fn setup<R: RngCore>(&mut self, rng: &mut R) -> (CRS<E>, Vec<Vec<E::ScalarField>>) {
        // Sample tau and compute its powers ==========================================================
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
        let scalar_size = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let g_table = FixedBase::get_window_table(scalar_size, window_size, g);
        let powers_of_g =
            FixedBase::msm::<E::G1>(scalar_size, window_size, &g_table, &powers_of_tau);

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

        let top_tau_table = FixedBase::get_window_table(scalar_size, window_size, g);
        let y =
            FixedBase::msm::<E::G1>(scalar_size, window_size, &top_tau_table, &top_tau);

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
        let lag_coeffs: Vec<E::ScalarField> = lagrange_coefficients(lag_domain, gamma);

        // save the long_term_secret L_B(gamma)
        self.long_term_secret = lag_coeffs[self.batch_size];

        // secret share the lagrange coefficients
        let mut lag_shares: Vec<Vec<E::ScalarField>> = Vec::new();
        let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(self.n).unwrap();
        for i in 0..self.batch_size {
            // secret share i-th coefficient
            let mut coeffs = vec![E::ScalarField::zero(); self.n];
            coeffs[0] = lag_coeffs[i];
            for j in 1..self.n/2 {
                coeffs[j] = E::ScalarField::rand(rng);
            }

            // use fft to compute shares
            let evals = share_domain.fft(&coeffs);
            lag_shares.push(evals);
        }
        
        let crs = CRS::<E> {
            powers_of_g,
            htau: h * tau,
            y,
        };

        (crs, transpose(lag_shares))
    }

    pub fn epoch_setup<R: RngCore>(&mut self, rng: &mut R) -> (E::G1, E::G2, E::G1, Vec<E::ScalarField>, Vec<E::ScalarField>) {
        // Generators
        let g = E::G1::generator();
        let h = E::G2::generator();

        // sample gtilde and htilde such that dlog_g(gtilde) = dlog_h(htilde)
        // publish gtilde and htilde
        let delta = E::ScalarField::rand(rng);
        let gtilde = g * delta;
        let htilde = h * delta;

        // secret share alpha*long_term_secret, r for a random alpha and r
        // publish com = g^alpha . gtilde^r
        let alpha = E::ScalarField::rand(rng);
        let r = E::ScalarField::rand(rng);
        let com = (g * alpha) + (gtilde * r);

        let mut alpha_coeffs = vec![E::ScalarField::zero(); self.n];
        let mut r_coeffs = vec![E::ScalarField::zero(); self.n];
        alpha_coeffs[0] = alpha * self.long_term_secret;
        r_coeffs[0] = r;
        for j in 1..self.n/2 {
            alpha_coeffs[j] = E::ScalarField::rand(rng);
            r_coeffs[j] = E::ScalarField::rand(rng);
        }

        // use fft to compute shares
        let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(self.n).unwrap();
        let alpha_shares = share_domain.fft(&alpha_coeffs);
        let r_shares = share_domain.fft(&r_coeffs);

        (gtilde, htilde, com, alpha_shares, r_shares)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    type E = Bls12_381;
    type G1 = <E as Pairing>::G1;
    type G2 = <E as Pairing>::G2;

    #[test]
    fn test_dealer() {
        let mut rng = ark_std::test_rng();
        let batch_size = 1 << 5;
        let n = 1 << 4;

        let mut dealer = Dealer::<E>::new(batch_size, n);
        let (crs, lag_shares) = dealer.setup(&mut rng);

        let (gtilde, htilde, _com, ashares, rshares) = dealer.epoch_setup(&mut rng);

        assert_eq!(E::pairing(gtilde, G2::generator()), E::pairing(G1::generator(), htilde));
        assert_eq!(lag_shares.len(), n);
        assert_eq!(lag_shares[0].len(), batch_size);
        assert_eq!(crs.powers_of_g.len(), batch_size + 1);
        assert_eq!(ashares.len(), n);
        assert_eq!(rshares.len(), n);
    }
}
