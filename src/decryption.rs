use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::*;
use ark_std::{end_timer, start_timer, Zero};

use crate::{encryption::Ciphertext, utils::hash_to_bytes};

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

    pub fn partial_decrypt(&self, ct: Vec<Ciphertext<E>>) -> E::ScalarField {
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
    ct: Vec<Ciphertext<E>>,
) {
    let batch_size = ct.len();
    let n = partial_decryptions.len();

    let share_domain = Radix2EvaluationDomain::<E::ScalarField>::new(n).unwrap();
    let tx_domain = Radix2EvaluationDomain::<E::ScalarField>::new(batch_size + 1).unwrap();

    let pof1 = share_domain.ifft(&partial_decryptions)[0];

    let mut pevals = vec![E::ScalarField::zero(); batch_size + 1];
    pevals[0] = pof1;

    //use FK22 to get all the KZG proofs quickly
}

// pub fn batch_open()
