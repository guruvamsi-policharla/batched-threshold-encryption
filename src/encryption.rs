use ark_ec::{pairing::Pairing, Group};
use ark_ff::Field;
use ark_serialize::*;
use ark_std::{rand::RngCore, UniformRand};
use merlin::Transcript;
use retry::{delay::NoDelay, retry};

use crate::utils::{hash_to_bytes, xor};

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct DLogProof<E: Pairing> {
    pub u: E::G1,
    pub z: E::ScalarField,
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
    pub fn verify(&self) -> bool {
        let g = E::G1::generator();

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
        self.pi.u.serialize_uncompressed(&mut u_bytes).unwrap();
        ts.append_message(&[6u8], &u_bytes);

        let mut c_bytes = Vec::new();
        ts.challenge_bytes(&[7u8], &mut c_bytes);
        let c = E::ScalarField::from_random_bytes(&c_bytes).unwrap();

        let lhs = g * self.pi.z;
        let rhs = self.pi.u + (self.gs * c);

        lhs == rhs
    }
}

pub fn encrypt<E: Pairing, R: RngCore>(
    msg: [u8; 32],
    x: E::ScalarField,
    com: E::G1,
    htilde: E::G2,
    pk: E::G2,
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
    let ct2 = (pk - (h * x)) * rho;
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

    let r = E::ScalarField::rand(rng);
    let u = g * r;
    let mut u_bytes = Vec::new();
    u.serialize_uncompressed(&mut u_bytes).unwrap();
    ts.append_message(&[6u8], &u_bytes);

    // Fiat-Shamir to get challenge
    let mut c_bytes = Vec::new();
    ts.challenge_bytes(&[7u8], &mut c_bytes);
    let c = E::ScalarField::from_random_bytes(&c_bytes).unwrap();

    // compute response
    let z = r + c * s;

    debug_assert_eq!(g * z, u + gs * c);

    let pi = DLogProof { u, z };

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
        let (_, htilde, com, _, _) = dealer.epoch_setup(&mut rng);

        let msg = [1u8; 32];
        let x = tx_domain.group_gen;

        let ct = encrypt::<Bls12_381, _>(msg, x, com, htilde, crs.pk, &mut rng);

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
        
        g.serialize_uncompressed(&mut g1_bytes).unwrap();
        h.serialize_uncompressed(&mut g2_bytes).unwrap();
        x.serialize_uncompressed(&mut fr_bytes).unwrap();

        println!("G1 len: {} bytes", g1_bytes.len());
        println!("G2 len: {} bytes", g2_bytes.len());
        println!("Fr len: {} bytes", fr_bytes.len());
    }
}
