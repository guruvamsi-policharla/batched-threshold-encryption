use ark_bls12_381::Bls12_381;
use ark_ec::{bls12::Bls12, pairing::Pairing};
use ark_poly::{Radix2EvaluationDomain, EvaluationDomain};
use batch_threshold::{dealer::Dealer, decryption::SecretKey, encryption::{encrypt, Ciphertext}};
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

type E = Bls12_381;
type Fr = <Bls12<ark_bls12_381::Config> as Pairing>::ScalarField;

//todo: use seeded randomness
fn bench_partial_decrypt(c: &mut Criterion) {
    let mut rng = ark_std::test_rng();
    
    let n = 1 << 4;
    let mut group = c.benchmark_group("partial_decrypt");
    for size in 2..=10 {
        let batch_size = 1 << size;   

        let mut dealer = Dealer::<E>::new(batch_size, n);
        let (crs, lag_shares) = dealer.setup(&mut rng);
        let (com, epoch_shares) = dealer.epoch_setup(&mut rng);

        let mut secret_key: Vec<SecretKey<E>> = Vec::new();
        for i in 0..n {
            secret_key.push(SecretKey::new(lag_shares[i].clone(), epoch_shares[i]));
        }

        let msg = [1u8; 32];

        let tx_domain = Radix2EvaluationDomain::<Fr>::new(batch_size).unwrap();
        
        // generate ciphertexts for all points in tx_domain
        let mut ct: Vec<Ciphertext<E>> = Vec::new();
        for x in tx_domain.elements() {
            ct.push(encrypt::<E, _>(msg, x, com, crs.pk, &mut rng));
        }

        // bench partial decryption
        group.bench_with_input(BenchmarkId::from_parameter(batch_size), &ct, |b, ct| {
            b.iter(|| secret_key[0].partial_decrypt(ct));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_partial_decrypt);
criterion_main!(benches);