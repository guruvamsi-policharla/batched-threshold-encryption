use ark_ec::{pairing::Pairing, Group};
use ark_ff::{batch_inversion, FftField, Field};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::*;
use ark_std::{end_timer, start_timer, One, Zero};

use crate::{
    dealer::CRS,
    encryption::Ciphertext,
    utils::{hash_to_bytes, xor, interpolate_almostgood},
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

pub fn public_reconstruction<E: Pairing>(
    partial_decryptions: Vec<E::ScalarField>,
    ct: &Vec<Ciphertext<E>>,
    crs: &CRS<E>,
) -> Vec<[u8; 32]> {
    let batch_size = ct.len();
    let n = partial_decryptions.len();

    let g = E::G1::generator();

    let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(n).unwrap();
    let top_domain = Radix2EvaluationDomain::<E::ScalarField>::new(2 * batch_size).unwrap();
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
    let mut v = f.clone();
    v.reverse();
    v.truncate(batch_size);
    v.resize(2 * batch_size, E::ScalarField::zero());
    let v = top_domain.fft(&v);

    // get all the roots of top_domain
    let top_domain_roots: Vec<E::ScalarField> = top_domain.elements().collect();

    // h = crs.powers_of_top_tau[i] ^ (v[i] . top_domain_roots[i])
    let mut h = vec![g; 2 * batch_size];
    for i in 0..2 * batch_size {
        h[i] = crs.powers_of_top_tau[i] * (v[i] * top_domain_roots[i]);
    }

    // inverse fft on h
    let mut h = top_domain.ifft(&h);
    h.resize(batch_size, E::G1::zero());

    // fft on h to get KZG proofs
    let pi = tx_domain.fft(&h);
    // pi.reverse();

    // now decrypt each of the ciphertexts as m = ct1 \xor H(e(pi, ct2))
    let mut m = vec![[0u8; 32]; batch_size];
    for i in 0..batch_size {
        let mask = E::pairing(pi[i], ct[i].ct2);
        let hmask = hash_to_bytes(mask);
        m[i] = xor(&ct[i].ct1, &hmask).as_slice().try_into().unwrap();
    }

    m
}
