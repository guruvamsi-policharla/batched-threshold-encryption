use ark_ec::{pairing::Pairing, Group};
use ark_ff::{FftField, Field};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::*;
use ark_std::Zero;

use crate::{
    dealer::CRS,
    encryption::Ciphertext,
    utils::{hash_to_bytes, interpolate_almostgood, open_all_values, xor},
};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct SecretKey<E: Pairing> {
    lag_share: Vec<E::ScalarField>,
    alpha_share: E::ScalarField,
    r_share: E::ScalarField,
}

impl<E: Pairing> SecretKey<E> {
    pub fn new(lag_share: Vec<E::ScalarField>, alpha_share: E::ScalarField, r_share: E::ScalarField) -> Self {
        SecretKey {
            lag_share,
            alpha_share,
            r_share,
        }
    }

    /// each party in the committee computes a partial decryption
    pub fn partial_decrypt(&self, ct: &Vec<Ciphertext<E>>) -> (E::ScalarField, E::G1) {

        let mut partial_decryption1 = self.alpha_share;
        let batch_size = self.lag_share.len();

        // Check that all ciphertexts are valid
        for i in 0..batch_size {
            assert!(ct[i].verify());
        }

        // compute partial decryption
        for i in 0..batch_size {
            let tg_bytes = hash_to_bytes(ct[i].gs);
            let peval = E::ScalarField::from_random_bytes(&tg_bytes).unwrap();
            partial_decryption1 += self.lag_share[i] * peval;
        }

        let partial_decryption2 =  E::G1::generator() * self.r_share;
        
        (partial_decryption1, partial_decryption2)
    }
}

/// decrypts all the ciphertexts in a batch
pub fn decrypt_all<E: Pairing>(
    partial_decryptions1: &Vec<E::ScalarField>,
    partial_decryptions2: &Vec<E::G1>,
    ct: &Vec<Ciphertext<E>>,
    crs: &CRS<E>,
) -> Vec<[u8; 32]> {
    let batch_size = ct.len();
    let n = partial_decryptions1.len();

    let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(n).unwrap();
    let tx_domain = Radix2EvaluationDomain::<E::ScalarField>::new(batch_size).unwrap();
    let gamma = E::ScalarField::GENERATOR;

    let fofgamma = share_domain.ifft(&partial_decryptions1)[0];
    let pi2 = share_domain.ifft(&partial_decryptions2)[0];

    // compute fevals by hashing gs of the ciphertexts to get fevals
    let mut fevals = vec![E::ScalarField::zero(); batch_size + 1];
    for i in 0..batch_size {
        let tg_bytes = hash_to_bytes(ct[i].gs);
        fevals[i] = E::ScalarField::from_random_bytes(&tg_bytes).unwrap();
    }
    fevals[batch_size] = fofgamma;

    // fevals are on an 'almost' nice domain. so we first interpolate quotient polynomial
    // where the evaluations are determined as q(x) = (f(x) - f(gamma))/(x-gamma)
    let f = interpolate_almostgood(&fevals, &tx_domain, fofgamma, gamma);

    // use FK22 to get all the KZG proofs in O(nlog n) time =======================
    let pi1 = open_all_values::<E>(&crs.y, &f, &tx_domain);

    // now decrypt each of the ciphertexts as m = ct1 \xor H(e(pi1, ct2).e(pi2,ct3))
    let mut m = vec![[0u8; 32]; batch_size];
    for i in 0..batch_size {
        let mask = E::pairing(pi1[i], ct[i].ct2) + E::pairing(pi2, ct[i].ct3);
        let hmask = hash_to_bytes(mask);
        m[i] = xor(&ct[i].ct1, &hmask).as_slice().try_into().unwrap();
    }

    m
}
