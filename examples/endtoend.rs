use ark_bls12_381::Bls12_381;
use ark_ec::pairing::Pairing;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};

use batch_threshold::{
    dealer::Dealer,
    decryption::{decrypt_all, SecretKey},
    encryption::{encrypt, Ciphertext},
};

type E = Bls12_381;
type Fr = <E as Pairing>::ScalarField;
type G1 = <E as Pairing>::G1;

fn main() {
    let mut rng = ark_std::test_rng();
    let batch_size = 1 << 5;
    let n = 1 << 4;

    let mut dealer = Dealer::<E>::new(batch_size, n);
    let (crs, lag_shares) = dealer.setup(&mut rng);
    let (gtilde, htilde, com, alpha_shares, r_shares) = dealer.epoch_setup(&mut rng);

    let mut secret_key: Vec<SecretKey<E>> = Vec::new();
    for i in 0..n {
        secret_key.push(SecretKey::new(lag_shares[i].clone(), alpha_shares[i], r_shares[i]));
    }

    let tx_domain = Radix2EvaluationDomain::<Fr>::new(batch_size).unwrap();

    let msg = [1u8; 32];

    // generate ciphertexts for all points in tx_domain
    let mut ct: Vec<Ciphertext<E>> = Vec::new();
    for x in tx_domain.elements() {
        ct.push(encrypt::<E, _>(msg, x, com, htilde, crs.htau, &mut rng));
    }

    // generate partial decryptions
    let mut partial_decryptions1: Vec<Fr> = Vec::new();
    let mut partial_decryptions2: Vec<G1> = Vec::new();
    for i in 0..n {
        let partial_decryption = secret_key[i].partial_decrypt(&ct, gtilde, htilde, &crs);
        partial_decryptions1.push(partial_decryption.0);
        partial_decryptions2.push(partial_decryption.1);
    }

    let messages = decrypt_all(&partial_decryptions1, &partial_decryptions2, &ct, &crs);
    for i in 0..batch_size {
        assert_eq!(msg, messages[i]);
    }
}
