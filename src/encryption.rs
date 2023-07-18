use ark_ec::{pairing::Pairing, Group};
use ark_std::{end_timer, rand::RngCore, start_timer, UniformRand};
use ark_serialize::*;
use ark_ff::Field;
use merlin::Transcript;

pub struct DLogProof<E: Pairing> {
    pub u: E::G1,
    pub z: E::ScalarField,
}

pub struct Ciphertext<E: Pairing> {
    pub ct1: [u8; 32],
    pub ct2: E::G1,
    pub gs: E::G1,
    pub x: E::ScalarField,
    pub pi: DLogProof<E>
}

pub fn encrypt<E:Pairing, R: RngCore>(msg: [u8; 32], x: E::ScalarField, com: E::G1, pk: E::G1, rng: &mut R) -> Ciphertext<E> {
    
    let enc_timer = start_timer!(|| "Encrypting"); 

    let g = E::G1::generator();
    let h = E::G2::generator();
    let rho = E::ScalarField::rand(rng);
    
    // hash element S to curve to get tg
    let s = E::ScalarField::rand(rng);
    let gs = g*s;
    let hgs = hash_to_bytes(gs);
    let tg = E::ScalarField::from_random_bytes(&hgs).unwrap();

    // compute mask
    let mask = E::pairing(com - (g*tg), h)*rho; //e(com/g^tg, h)^rho
    let hmask = hash_to_bytes(mask);

    // xor msg and hmask
    let mut ct1 = [0u8; 32];
    for i in 0..32 {
        ct1[i] = msg[i] ^ hmask[i];
    }

    let ct2 = (pk - g*x)*rho;

    // Prove knowledge of discrete log of S with ct1, ct2, S, x as tags
    let mut ts: Transcript = Transcript::new(&[0u8]);
    ts.append_message(&[1u8], &ct1);
    
    let mut ct2_bytes = Vec::new(); ct2.serialize_uncompressed(&mut ct2_bytes).unwrap();
    ts.append_message(&[2u8], &ct2_bytes);

    let mut gs_bytes = Vec::new(); gs.serialize_uncompressed(&mut gs_bytes).unwrap();
    ts.append_message(&[3u8], &gs_bytes);

    let mut x_bytes = Vec::new(); x.serialize_uncompressed(&mut x_bytes).unwrap();
    ts.append_message(&[4u8], &x_bytes);

    let r = E::ScalarField::rand(rng);
    let u = g*r;
    let mut u_bytes = Vec::new(); u.serialize_uncompressed(&mut u_bytes).unwrap();
    ts.append_message(&[5u8], &u_bytes);

    // Fiat-Shamir to get challenge
    let mut c_bytes = Vec::new();
    ts.challenge_bytes(&[6u8], &mut c_bytes);
    let c = E::ScalarField::from_random_bytes(&c_bytes).unwrap();

    // compute response
    let z = r + c*s;

    debug_assert_eq!(g*z, u + gs*c);

    let pi = DLogProof {
        u,
        z
    };
    
    end_timer!(enc_timer);

    Ciphertext {
        ct1,
        ct2,
        gs,
        x,
        pi
    }
}

fn hash_to_bytes<T: CanonicalSerialize>(inp: T) -> [u8; 32] {
    let mut bytes = Vec::new();
    inp.serialize_uncompressed(&mut bytes).unwrap();
    let hash = blake3::hash(bytes.as_slice());
    let hash_bytes = hash.as_bytes();
    *hash_bytes
}