use ark_bls12_381::Bls12_381;
use ark_ec::{pairing::Pairing, bls12::Bls12};
use ark_poly::{Radix2EvaluationDomain, EvaluationDomain};
use batch_threshold::{dealer::Dealer, encryption::{encrypt, Ciphertext}, decryption::{SecretKey, public_reconstruction}};

type E = Bls12_381;
type Fr = <Bls12<ark_bls12_381::Config> as Pairing>::ScalarField;
// type G1 = <Bls12<ark_bls12_381::Config> as Pairing>::G1;
// type G2 = <Bls12<ark_bls12_381::Config> as Pairing>::G2;

fn main() {
    let mut rng = ark_std::test_rng();
    let batch_size = 1 << 5;
    let n = 1 << 4;

    // let g = G1::generator();
    // let h = G2::generator();

    let mut dealer = Dealer::<E>::new(batch_size, n);
    let (crs, lag_shares) = dealer.setup(&mut rng);

    let (com, evals) = dealer.epoch_setup(&mut rng);

    let mut secret_key: Vec<SecretKey<E>> = Vec::new();
    for i in 0..n {
        secret_key.push(SecretKey::new(lag_shares[i].clone(), evals[i]));
    }

    let tx_domain = Radix2EvaluationDomain::<Fr>::new(batch_size).unwrap();

    let msg = [1u8; 32];

    // generate ciphertexts for all points in tx_domain
    let mut ct:Vec<Ciphertext<E>> = Vec::new();
    for x in tx_domain.elements() {
        ct.push(encrypt::<E, _>(msg, x, com, crs.pk, &mut rng));
    }

    // generate partial decryptions
    let mut partial_decryptions: Vec<Fr> = Vec::new();
    for i in 0..n {
        partial_decryptions.push(secret_key[i].partial_decrypt(&ct));
    }

    // interpolate the polynomial
    let messages = public_reconstruction(partial_decryptions, &ct, &crs);

    println!("messages: {:?}", messages);
}