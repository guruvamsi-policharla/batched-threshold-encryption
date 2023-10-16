use ark_bls12_381::Bls12_381;
use ark_ec::pairing::Pairing;
use ark_poly::{Radix2EvaluationDomain, EvaluationDomain};
use batch_threshold::{dealer::Dealer, decryption::SecretKey, encryption::{encrypt, Ciphertext}};
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

type E = Bls12_381;
type Fr = <E as Pairing>::ScalarField;

//todo: use seeded randomness
fn bench_partial_decrypt(c: &mut Criterion) {
    let mut rng = ark_std::test_rng();
    
    let n = 1 << 4;
    let mut group = c.benchmark_group("partial_decrypt");
    group.sample_size(20);
    
    for size in 2..=10 {
        let batch_size = 1 << size;   

        let mut dealer = Dealer::<E>::new(batch_size, n);
        let (crs, lag_shares) = dealer.setup(&mut rng);
        let (gtilde, htilde, com, alpha_shares, r_shares) = dealer.epoch_setup(&mut rng);

        let mut secret_key: Vec<SecretKey<E>> = Vec::new();
        for i in 0..n {
            secret_key.push(SecretKey::new(lag_shares[i].clone(), alpha_shares[i], r_shares[i]));
        }

        let msg = [1u8; 32];

        let tx_domain = Radix2EvaluationDomain::<Fr>::new(batch_size).unwrap();
        
        // generate ciphertexts for all points in tx_domain
        let mut ct: Vec<Ciphertext<E>> = Vec::new();
        for x in tx_domain.elements() {
            ct.push(encrypt::<E, _>(msg, x, com, htilde, crs.htau, &mut rng));
        }

        // bench partial decryption
        group.bench_with_input(BenchmarkId::from_parameter(batch_size), &ct, |b, ct| {
            b.iter(|| secret_key[0].partial_decrypt(&ct, gtilde, htilde, &crs));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_partial_decrypt);
criterion_main!(benches);