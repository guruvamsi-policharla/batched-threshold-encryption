use ark_ff::{batch_inversion, FftField};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use ark_serialize::CanonicalSerialize;

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
    let mut den: Vec<F> = vec![F::one(); domain.size()];
    let omega = domain.group_gen;
    for i in 1..domain.size() {
        den[i] = den[i - 1] * omega;
    }
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

#[cfg(test)]
mod tests {
    use ark_ec::{bls12::Bls12, pairing::Pairing};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::{UniformRand, Zero};

    use super::*;
    type Fr = <Bls12<ark_bls12_381::Config> as Pairing>::ScalarField;

    #[test]
    fn lagrange_test() {
        let mut rng = ark_std::test_rng();

        let domain_size = 1 << 5;
        let domain = Radix2EvaluationDomain::<Fr>::new(domain_size).unwrap();

        // almost_domain is {1, omega, ..., omega^(domain_size-1), tau}
        let tau = Fr::rand(&mut rng);
        let mut almost_domain: Vec<Fr> = domain.elements().collect();
        almost_domain.resize(domain_size+1, tau);

        let gamma = Fr::GENERATOR;
        let lag_coeffs = lagrange_coefficients(almost_domain, gamma);

        // sample a random vector of field elements as evals
        let evals: Vec<Fr> = (0..domain_size).map(|_| Fr::rand(&mut rng)).collect();
        let foftau = Fr::rand(&mut rng);

        let mut fofgamma = Fr::zero();
        for i in 0..domain_size {
            fofgamma += lag_coeffs[i] * evals[i];
        }
        fofgamma += lag_coeffs[domain_size]*foftau;

        let f = DensePolynomial::from_coefficients_vec(interpolate_almostgood(
            &evals,
            &domain,
            foftau,
            tau,
        ));

        assert_eq!(f.evaluate(&gamma), fofgamma);
    }
}
