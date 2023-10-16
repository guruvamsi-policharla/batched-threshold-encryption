use ark_bls12_381::Bls12_381;
use ark_ec::pairing::Pairing;
use batch_threshold::{dealer::Dealer, decryption::SecretKey, encryption::encrypt};
use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use ark_std::One;

type E = Bls12_381;
type Fr = <E as Pairing>::ScalarField;

//todo: use seeded randomness
fn bench_encrypt(c: &mut Criterion) {
    let mut rng = ark_std::test_rng();
    
    let n = 1 << 4;
    let mut group = c.benchmark_group("encrypt");
    for size in 2..6 {
        // timing doesn't change since ecryption is independent of batch_size
        // done as a sanity check

        let batch_size = 1 << size;   

        let mut dealer = Dealer::<E>::new(batch_size, n);
        let (crs, lag_shares) = dealer.setup(&mut rng);
        let (_gtilde, htilde, com, alpha_shares, r_shares) = dealer.epoch_setup(&mut rng);

        let mut secret_key: Vec<SecretKey<E>> = Vec::new();
        for i in 0..n {
            secret_key.push(SecretKey::new(lag_shares[i].clone(), alpha_shares[i], r_shares[i]));
        }

        let msg = [1u8; 32];

        group.bench_with_input(BenchmarkId::from_parameter(batch_size), &(msg, com, crs.htau), |b, &inp| {
            b.iter(|| encrypt::<E, _>(inp.0, Fr::one(), inp.1, htilde, inp.2, &mut rng));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_encrypt);
criterion_main!(benches);