use ark_ec::{pairing::Pairing, Group};
use ark_ff::Field;
use ark_serialize::*;
use ark_std::{rand::RngCore, UniformRand};
use merlin::Transcript;
use retry::{delay::NoDelay, retry};

use crate::{utils::{hash_to_bytes, xor}, dealer::CRS};

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct DLogProof<E: Pairing> {
    pub c: E::ScalarField, //challenge
    pub z1: E::ScalarField, //opening for g^s
    pub z2: E::ScalarField, //opening for htilde^{rho}
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct Ciphertext<E: Pairing> {
    pub ct1: [u8; 32],
    pub ct2: E::G2,
    pub ct3: E::G2,
    pub gs: E::G1,
    pub x: E::ScalarField,
    pub pi: DLogProof<E>,
}

impl<E: Pairing> Ciphertext<E> {
    /// panicks if ciphertext does not verify
    pub fn verify(&self, gtilde: E::G1, htilde: E::G2, crs: &CRS<E>) {
        let g = E::G1::generator();

        let u = g * self.pi.z1 - self.gs * self.pi.c;
        let v = htilde * self.pi.z2 - self.ct3 * self.pi.c;
        
        let mut ts: Transcript = Transcript::new(&[0u8]);
        ts.append_message(&[1u8], &self.ct1);

        let mut ct2_bytes = Vec::new();
        self.ct2.serialize_uncompressed(&mut ct2_bytes).unwrap();
        ts.append_message(&[2u8], &ct2_bytes);

        let mut ct3_bytes = Vec::new();
        self.ct3.serialize_uncompressed(&mut ct3_bytes).unwrap();
        ts.append_message(&[3u8], &ct3_bytes);

        let mut gs_bytes = Vec::new();
        self.gs.serialize_uncompressed(&mut gs_bytes).unwrap();
        ts.append_message(&[4u8], &gs_bytes);

        let mut x_bytes = Vec::new();
        self.x.serialize_uncompressed(&mut x_bytes).unwrap();
        ts.append_message(&[5u8], &x_bytes);

        let mut u_bytes = Vec::new();
        u.serialize_uncompressed(&mut u_bytes).unwrap();
        ts.append_message(&[6u8], &u_bytes);
        
        let mut v_bytes = Vec::new();
        v.serialize_uncompressed(&mut v_bytes).unwrap();
        ts.append_message(&[7u8], &v_bytes);
        
        // Fiat-Shamir to get challenge
        let mut c_bytes = [0u8;31];
        ts.challenge_bytes(&[8u8], &mut c_bytes);
        let c = E::ScalarField::from_random_bytes(&c_bytes).unwrap();

        // assert that the recomputed challenge matches
        assert_eq!(self.pi.c, c);

        // pairing check
        let gtaux = crs.powers_of_g[1] - (g*self.x);
        assert_eq!(
            E::pairing(gtilde, self.ct2),
            E::pairing(gtaux, self.ct3)
        );        
    }
}

pub fn encrypt<E: Pairing, R: RngCore>(
    msg: [u8; 32],
    x: E::ScalarField,
    com: E::G1,
    htilde: E::G2,
    htau: E::G2,
    rng: &mut R,
) -> Ciphertext<E> {
    let g = E::G1::generator();
    let h = E::G2::generator();
    let rho = E::ScalarField::rand(rng);

    // hash element S to curve to get tg
    // retry if bytes cannot be converted to a field element
    let result = retry(NoDelay, || {
        let s = E::ScalarField::rand(rng);
        let gs = g * s;
        let hgs = hash_to_bytes(gs);
        let tg_option = E::ScalarField::from_random_bytes(&hgs);

        match tg_option {
            Some(tg) => Ok((s, gs, tg)),
            None => {
                #[cfg(debug_assertions)]
                {
                    dbg!("Failed to hash to field element, retrying...");
                }
                Err(())
            }
        }
    });

    let (s, gs, tg) = result.unwrap();

    // compute mask
    let mask = E::pairing(com - (g * tg), h) * rho; //e(com/g^tg, h)^rho
    let hmask = hash_to_bytes(mask);

    // xor msg and hmask
    let ct1: [u8; 32] = xor(&msg, &hmask).as_slice().try_into().unwrap();
    let ct2 = (htau - (h * x)) * rho;
    let ct3 = htilde * rho;

    // Prove knowledge of discrete log of S with ct1, ct2, S, x as tags
    let mut ts: Transcript = Transcript::new(&[0u8]);
    ts.append_message(&[1u8], &ct1);

    let mut ct2_bytes = Vec::new();
    ct2.serialize_uncompressed(&mut ct2_bytes).unwrap();
    ts.append_message(&[2u8], &ct2_bytes);

    let mut ct3_bytes = Vec::new();
    ct3.serialize_uncompressed(&mut ct3_bytes).unwrap();
    ts.append_message(&[3u8], &ct3_bytes);

    let mut gs_bytes = Vec::new();
    gs.serialize_uncompressed(&mut gs_bytes).unwrap();
    ts.append_message(&[4u8], &gs_bytes);

    let mut x_bytes = Vec::new();
    x.serialize_uncompressed(&mut x_bytes).unwrap();
    ts.append_message(&[5u8], &x_bytes);

    let r1 = E::ScalarField::rand(rng);
    let r2 = E::ScalarField::rand(rng);
    let u = g * r1;
    let v = htilde * r2;
    
    let mut u_bytes = Vec::new();
    u.serialize_uncompressed(&mut u_bytes).unwrap();
    ts.append_message(&[6u8], &u_bytes);

    let mut v_bytes = Vec::new();
    v.serialize_uncompressed(&mut v_bytes).unwrap();
    ts.append_message(&[7u8], &v_bytes);

    // Fiat-Shamir to get challenge
    // note we sample a 31-byte field element to avoid rejection sampling. 
    // this is secure because the challenge only affects soundness, not zk.
    // the probability of a bad challenge is upperbounded by 1/2^{31*8}
    let mut c_bytes = [0u8;31];
    ts.challenge_bytes(&[8u8], &mut c_bytes);
    let c = E::ScalarField::from_random_bytes(&c_bytes).unwrap();

    // compute response
    let z1 = r1 + c * s;
    let z2 = r2 + c * rho;

    debug_assert_eq!(g * z1, u + gs * c);
    debug_assert_eq!(htilde * z2, v + ct3 * c);

    let pi = DLogProof { c, z1, z2 };

    Ciphertext {
        ct1,
        ct2,
        ct3,
        gs,
        x,
        pi,
    }
}

#[cfg(test)]
mod tests {
    use crate::dealer::Dealer;

    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_ec::bls12::Bls12;
    use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};

    type E = Bls12_381;
    type Fr = <Bls12<ark_bls12_381::Config> as Pairing>::ScalarField;
    type G1 = <Bls12<ark_bls12_381::Config> as Pairing>::G1;
    type G2 = <Bls12<ark_bls12_381::Config> as Pairing>::G2;

    #[test]
    fn test_encryption() {
        let mut rng = ark_std::test_rng();

        let batch_size = 1 << 5;
        let tx_domain = Radix2EvaluationDomain::<Fr>::new(batch_size).unwrap();

        let mut dealer = Dealer::<E>::new(batch_size, 1 << 4);
        let (crs, _) = dealer.setup(&mut rng);
        let (gtilde, htilde, com, _, _) = dealer.epoch_setup(&mut rng);

        let msg = [1u8; 32];
        let x = tx_domain.group_gen;

        let ct = encrypt::<Bls12_381, _>(msg, x, com, htilde, crs.htau, &mut rng);

        let mut ct_bytes = Vec::new();
        ct.serialize_compressed(&mut ct_bytes).unwrap();
        println!("Compressed ciphertext: {} bytes", ct_bytes.len());

        let mut ct_bytes = Vec::new();
        ct.serialize_uncompressed(&mut ct_bytes).unwrap();
        println!("Uncompressed ciphertext: {} bytes", ct_bytes.len());

        let mut g1_bytes = Vec::new();
        let mut g2_bytes = Vec::new();
        let mut fr_bytes = Vec::new();
        
        let g = G1::generator();
        let h = G2::generator();
        let x = tx_domain.group_gen;
        
        g.serialize_compressed(&mut g1_bytes).unwrap();
        h.serialize_compressed(&mut g2_bytes).unwrap();
        x.serialize_compressed(&mut fr_bytes).unwrap();

        println!("G1 len: {} bytes", g1_bytes.len());
        println!("G2 len: {} bytes", g2_bytes.len());
        println!("Fr len: {} bytes", fr_bytes.len());

        ct.verify(gtilde, htilde, &crs);
    }
}
