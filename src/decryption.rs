use ark_ec::pairing::Pairing;
use ark_ff::{FftField, Field};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::*;
use ark_std::{end_timer, start_timer, Zero};

use crate::{
    dealer::CRS,
    encryption::Ciphertext,
    utils::{hash_to_bytes, interpolate_almostgood, open_all_values, xor},
};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct SecretKey<E: Pairing> {
    lag_share: Vec<E::ScalarField>,
    epoch_share: E::ScalarField,
}

impl<E: Pairing> SecretKey<E> {
    pub fn new(lag_share: Vec<E::ScalarField>, epoch_share: E::ScalarField) -> Self {
        SecretKey {
            lag_share,
            epoch_share,
        }
    }

    /// each party in the committee computes a partial decryption
    pub fn partial_decrypt(&self, ct: &Vec<Ciphertext<E>>) -> E::ScalarField {
        let partial_dec_timer = start_timer!(|| "Partial decryption");

        let mut partial_decryption = self.epoch_share;
        let batch_size = self.lag_share.len();

        // Check that all ciphertexts are valid
        for i in 0..batch_size {
            assert!(ct[i].verify());
        }

        // compute partial decryption
        for i in 0..batch_size {
            let tg_bytes = hash_to_bytes(ct[i].gs);
            let peval = E::ScalarField::from_random_bytes(&tg_bytes).unwrap();
            partial_decryption += self.lag_share[i] * peval;
        }

        end_timer!(partial_dec_timer);

        partial_decryption
    }
}

/// decrypts all the ciphertexts in a batch
pub fn decrypt_all<E: Pairing>(
    partial_decryptions: Vec<E::ScalarField>,
    ct: &Vec<Ciphertext<E>>,
    crs: &CRS<E>,
) -> Vec<[u8; 32]> {
    let batch_size = ct.len();
    let n = partial_decryptions.len();

    let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(n).unwrap();
    let tx_domain = Radix2EvaluationDomain::<E::ScalarField>::new(batch_size).unwrap();
    let gamma = E::ScalarField::GENERATOR;

    let fofgamma = share_domain.ifft(&partial_decryptions)[0];

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
    let pi = open_all_values::<E>(&crs.y, &f, &tx_domain);

    // now decrypt each of the ciphertexts as m = ct1 \xor H(e(pi, ct2))
    let mut m = vec![[0u8; 32]; batch_size];
    for i in 0..batch_size {
        let mask = E::pairing(pi[i], ct[i].ct2);
        let hmask = hash_to_bytes(mask);
        m[i] = xor(&ct[i].ct1, &hmask).as_slice().try_into().unwrap();
    }

    m
}
