use ark_bls12_381::Bls12_381;
use ark_ec::pairing::Pairing;
use ark_poly::{Radix2EvaluationDomain, EvaluationDomain};
use batch_threshold::{dealer::Dealer, decryption::{SecretKey, decrypt_all}, encryption::{encrypt, Ciphertext}};
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

type E = Bls12_381;
type Fr = <E as Pairing>::ScalarField;
type G1 = <E as Pairing>::G1;

//todo: use seeded randomness
fn bench_decrypt_all(c: &mut Criterion) {
    let mut rng = ark_std::test_rng();
    
    let n = 1 << 4;
    let mut group = c.benchmark_group("decrypt_all");
    for size in 2..=10 {
        let batch_size = 1 << size;   

        let mut dealer = Dealer::<E>::new(batch_size, n);
        let (crs, lag_shares) = dealer.setup(&mut rng);
        let (_gtilde, htilde, com, alpha_shares, r_shares) = dealer.epoch_setup(&mut rng);

        let mut secret_key: Vec<SecretKey<E>> = Vec::new();
        for i in 0..n {
            secret_key.push(SecretKey::new(lag_shares[i].clone(), alpha_shares[i], r_shares[i]));
        }

        let msg = [1u8; 32];

        let tx_domain = Radix2EvaluationDomain::<Fr>::new(batch_size).unwrap();
        
        // generate ciphertexts for all points in tx_domain
        let mut ct: Vec<Ciphertext<E>> = Vec::new();
        for x in tx_domain.elements() {
            ct.push(encrypt::<E, _>(msg, x, com, htilde, crs.pk, &mut rng));
        }

        // generate partial decryptions
        let mut partial_decryptions1: Vec<Fr> = Vec::new();
        let mut partial_decryptions2: Vec<G1> = Vec::new();
        for i in 0..n {
            let partial_decryption = secret_key[i].partial_decrypt(&ct);
            partial_decryptions1.push(partial_decryption.0);
            partial_decryptions2.push(partial_decryption.1);
        }
        
        // bench full decryption
        group.bench_with_input(BenchmarkId::from_parameter(batch_size), &(partial_decryptions1, partial_decryptions2, ct, crs), |b, inp| {
            b.iter(|| decrypt_all(&inp.0, &inp.1, &inp.2, &inp.3));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_decrypt_all);
criterion_main!(benches);