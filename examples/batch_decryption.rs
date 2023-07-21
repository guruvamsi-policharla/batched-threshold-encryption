use ark_bls12_381::Bls12_381;
use ark_ec::{bls12::Bls12, pairing::Pairing, Group};
use ark_ff::{fields::Field, FftField};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial,
    Radix2EvaluationDomain,
};
use ark_std::Zero;
use batch_threshold::{
    dealer::Dealer,
    decryption::{public_reconstruction, SecretKey},
    encryption::{encrypt, Ciphertext},
    utils::{compute_opening_proof, hash_to_bytes, interpolate_almostgood, xor},
};

type E = Bls12_381;
type Fr = <Bls12<ark_bls12_381::Config> as Pairing>::ScalarField;
type G1 = <Bls12<ark_bls12_381::Config> as Pairing>::G1;
type G2 = <Bls12<ark_bls12_381::Config> as Pairing>::G2;

fn main() {
    let mut rng = ark_std::test_rng();
    let batch_size = 1 << 5;
    let n = 1 << 4;

    let g = G1::generator();
    let h = G2::generator();

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
    let mut ct: Vec<Ciphertext<E>> = Vec::new();
    for x in tx_domain.elements() {
        ct.push(encrypt::<E, _>(msg, x, com, crs.pk, &mut rng));
    }

    // generate partial decryptions
    let mut partial_decryptions: Vec<Fr> = Vec::new();
    for i in 0..n {
        partial_decryptions.push(secret_key[i].partial_decrypt(&ct));
    }

    // interpolate the polynomial with gamma
    let share_domain = Radix2EvaluationDomain::<Fr>::new(n).unwrap();
    let gamma = Fr::GENERATOR;
    let fofgamma = share_domain.ifft(&partial_decryptions)[0];

    // compute fevals by hashing gs of the ciphertexts to get fevals
    let mut fevals = vec![Fr::zero(); batch_size + 1];
    for i in 0..batch_size {
        let tg_bytes = hash_to_bytes(ct[i].gs);
        fevals[i] = Fr::from_random_bytes(&tg_bytes).unwrap();
    }
    fevals[batch_size] = fofgamma;

    let f = interpolate_almostgood(&fevals, &tx_domain, fofgamma, gamma);
    
    let fpoly = DensePolynomial::from_coefficients_vec(f.clone());

    // check that interpolation was carried out correctly
    for i in 0..batch_size {
        debug_assert_eq!(fpoly.evaluate(&tx_domain.element(i)), fevals[i]);
    }
    debug_assert_eq!(fpoly.evaluate(&gamma), fofgamma);

    // generate opening proofs
    // check that the kzg proof verifies
    for i in 0..batch_size {
        let pi = compute_opening_proof::<E>(&crs, &fpoly, &tx_domain.element(i));
        let lhs = E::pairing(com - (g * fpoly.evaluate(&tx_domain.element(i))), h);
        let rhs = E::pairing(pi, crs.pk - (h * tx_domain.element(i)));
        assert_eq!(lhs, rhs);
        // now decrypt each of the ciphertexts as m = ct1 \xor H(e(pi, ct2))
        let mask = E::pairing(pi, ct[i].ct2);
        let hmask = hash_to_bytes(mask);
        let m = xor(&ct[i].ct1, &hmask);

        // println!("m: {:?}", m);
    }
    
    let messages = public_reconstruction(partial_decryptions, &ct, &crs);
    // println!("messages: {:?}", messages);
}
