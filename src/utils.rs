use ark_ec::{pairing::Pairing, CurveGroup, VariableBaseMSM};
use ark_ff::{batch_inversion, FftField, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial,
    Radix2EvaluationDomain,
};
use ark_serialize::CanonicalSerialize;
use ark_std::{One, Zero};
use std::ops::Div;

use crate::dealer::CRS;

pub fn transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>> {
    assert!(!v.is_empty());
    let len = v[0].len();
    let mut iters: Vec<_> = v.into_iter().map(|n| n.into_iter()).collect();
    (0..len)
        .map(|_| {
            iters
                .iter_mut()
                .map(|n| n.next().unwrap())
                .collect::<Vec<T>>()
        })
        .collect()
}

pub fn hash_to_bytes<T: CanonicalSerialize>(inp: T) -> [u8; 32] {
    let mut bytes = Vec::new();
    inp.serialize_uncompressed(&mut bytes).unwrap();
    let hash = blake3::hash(bytes.as_slice());
    let hash_bytes = hash.as_bytes();
    *hash_bytes
}

pub fn xor(a: &[u8], b: &[u8]) -> Vec<u8> {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(x, y)| x ^ y).collect()
}

/// interpolates a polynomial on a domain that is *almost* good
/// say there exists a smooth domain of size d, and the degree of the polynomial is d
/// then there is one point outside the nice domain. to interpolate we use the fact
/// that the quotient polynomial is degree d-1 and interpolate it on the nice domain
/// and then recvoer the original polynomial
/// q(x) = (f(x) - f(gamma)) / (x-gamma), where gamma is the bad point
pub fn interpolate_almostgood<F: FftField>(
    fevals: &Vec<F>,
    domain: &Radix2EvaluationDomain<F>,
    fofgamma: F,
    gamma: F,
) -> Vec<F> {
    let mut qevals = vec![F::zero(); domain.size()];
    let mut den: Vec<F> = domain.elements().collect();
    den.iter_mut().for_each(|x| *x -= gamma);

    // batch invert den
    batch_inversion(&mut den);

    for i in 0..domain.size() {
        qevals[i] = (fevals[i] - fofgamma) * den[i];
    }
    let q = domain.ifft(&qevals);

    let mut f = vec![F::zero(); domain.size() + 1];
    f[0] = fofgamma - gamma * q[0];
    for i in 1..domain.size() {
        f[i] = q[i - 1] - gamma * q[i];
    }
    f[domain.size()] = q[domain.size() - 1];

    f
}

/// takes an input a domain and an evaluation point and generates lagrange coefficients for the same
pub fn lagrange_coefficients<F: FftField>(domain: Vec<F>, x: F) -> Vec<F> {
    let mut lag_coeffs: Vec<F> = vec![F::one(); domain.len()];
    for i in 0..domain.len() {
        let mut num = F::one();
        let mut den = F::one();
        for j in 0..domain.len() {
            if i != j {
                num *= x - domain[j];
                den *= domain[i] - domain[j];
            }
        }
        lag_coeffs[i] = num / den;
    }

    lag_coeffs
}

/// compute KZG opening proof
pub fn compute_opening_proof<E: Pairing>(
    crs: &CRS<E>,
    polynomial: &DensePolynomial<E::ScalarField>,
    point: &E::ScalarField,
) -> E::G1 {
    let eval = polynomial.evaluate(point);
    let eval_as_poly = DensePolynomial::from_coefficients_vec(vec![eval]);
    let numerator = polynomial - &eval_as_poly;
    let divisor = DensePolynomial::from_coefficients_vec(vec![
        E::ScalarField::zero() - point,
        E::ScalarField::one(),
    ]);
    let witness_polynomial = numerator.div(&divisor);

    commit_g1::<E>(&crs.powers_of_g, &witness_polynomial)
}

pub fn commit_g1<E: Pairing>(srs: &[E::G1], polynomial: &DensePolynomial<E::ScalarField>) -> E::G1 {
    if srs.len() - 1 < polynomial.degree() {
        panic!(
            "SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
            polynomial.degree(),
            srs.len()
        );
    }

    let plain_coeffs = convert_to_bigints(&polynomial.coeffs());
    let affine_srs = srs
        .iter()
        .map(|g| g.into_affine())
        .collect::<Vec<E::G1Affine>>();
    <E::G1 as VariableBaseMSM>::msm_bigint(&affine_srs, &plain_coeffs)
}

fn convert_to_bigints<F: PrimeField>(p: &[F]) -> Vec<F::BigInt> {
    let coeffs = ark_std::cfg_iter!(p)
        .map(|s| s.into_bigint())
        .collect::<Vec<_>>();
    coeffs
}

/// Computes all the openings of a KZG commitment in O(n log n) time
/// See https://github.com/khovratovich/Kate/blob/master/Kate_amortized.pdf
/// eprint version has a bug and hasn't been updated
pub fn open_all_values<E: Pairing>(
    powers_of_g: &Vec<E::G1>,
    f: &Vec<E::ScalarField>,
    domain: &Radix2EvaluationDomain<E::ScalarField>,
) -> Vec<E::G1> {
    let top_domain = Radix2EvaluationDomain::<E::ScalarField>::new(2 * domain.size()).unwrap();

    // use FK22 to get all the KZG proofs in O(nlog n) time =======================
    // todo: preprocess below
    let mut y = powers_of_g.clone();
    y.truncate(domain.size());
    y.reverse();
    y.resize(2 * domain.size(), E::G1::zero());
    let y = top_domain.fft(&y);

    // if f = {f0 ,f1, ..., fd}
    // v = {fd, (d-1 0s), fd, f, f1, ..., fd-2}
    let f = f[1..f.len()].to_vec();
    let mut v = vec![f[f.len() - 1]];
    v.resize(domain.size(), E::ScalarField::zero());
    v.push(f[f.len() - 1]);
    for &e in f.iter().take(f.len() - 1) {
        v.push(e);
    }
    assert_eq!(v.len(), 2 * domain.size());
    let v = top_domain.fft(&v);

    // h = y \odot v
    let mut h = vec![E::G1::zero(); 2 * domain.size()];
    for i in 0..2 * domain.size() {
        h[i] = y[i] * (v[i]);
    }

    // inverse fft on h
    let mut h = top_domain.ifft(&h);

    h.truncate(domain.size());

    // fft on h to get KZG proofs
    let pi = domain.fft(&h);

    pi
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_ec::{bls12::Bls12, pairing::Pairing, Group};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::{UniformRand, Zero};

    use crate::dealer::Dealer;

    use super::*;
    type Fr = <Bls12<ark_bls12_381::Config> as Pairing>::ScalarField;
    type G1 = <Bls12<ark_bls12_381::Config> as Pairing>::G1;
    type G2 = <Bls12<ark_bls12_381::Config> as Pairing>::G2;
    type E = Bls12_381;

    #[test]
    fn open_all_test() {
        let mut rng = ark_std::test_rng();

        let domain_size = 1 << 5;
        let domain = Radix2EvaluationDomain::<Fr>::new(domain_size).unwrap();

        let mut dealer = Dealer::<E>::new(domain_size, 1 << 5);
        let (crs, _) = dealer.setup(&mut rng);

        let mut f = vec![Fr::zero(); domain_size + 1];
        for i in 0..domain_size {
            f[i] = Fr::rand(&mut rng);
        }
        let fpoly = DensePolynomial::from_coefficients_vec(f.clone());
        let com = commit_g1::<E>(&crs.powers_of_g, &fpoly);
        let pi = open_all_values::<E>(&crs.powers_of_g, &f, &domain);

        // verify the kzg proof
        let g = G1::generator();
        let h = G2::generator();

        for i in 0..domain_size {
            let lhs = E::pairing(com - (g * fpoly.evaluate(&domain.element(i))), h);
            let rhs = E::pairing(pi[i], crs.pk - (h * domain.element(i)));
            assert_eq!(lhs, rhs);
        }
    }

    #[test]
    fn lagrange_test() {
        let mut rng = ark_std::test_rng();

        let domain_size = 1 << 5;
        let domain = Radix2EvaluationDomain::<Fr>::new(domain_size).unwrap();

        // almost_domain is {1, omega, ..., omega^(domain_size-1), tau}
        let tau = Fr::rand(&mut rng);
        let mut almost_domain: Vec<Fr> = domain.elements().collect();
        almost_domain.resize(domain_size + 1, tau);

        let gamma = Fr::GENERATOR;
        let lag_coeffs = lagrange_coefficients(almost_domain, gamma);

        // sample a random vector of field elements as evals
        let evals: Vec<Fr> = (0..domain_size).map(|_| Fr::rand(&mut rng)).collect();
        let foftau = Fr::rand(&mut rng);

        let mut fofgamma = Fr::zero();
        for i in 0..domain_size {
            fofgamma += lag_coeffs[i] * evals[i];
        }
        fofgamma += lag_coeffs[domain_size] * foftau;

        let f = DensePolynomial::from_coefficients_vec(interpolate_almostgood(
            &evals, &domain, foftau, tau,
        ));

        assert_eq!(f.evaluate(&gamma), fofgamma);
    }
}
