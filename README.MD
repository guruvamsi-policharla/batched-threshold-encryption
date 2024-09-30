# ‼️UPDATE: For virtually all applications, the [improved batched-threshold encryption scheme](https://github.com/guruvamsi-policharla/batched-threshold-pp) should be preferred.
# Batched Threshold Encryption [ePrint:2024/669](https://eprint.iacr.org/2024/669)

Rust implementation of the batched-threshold encryption introduced in [ePrint:2024/669](https://eprint.iacr.org/2024/669).

Use ```cargo bench``` to benchmark `encrypt`, `partial_decrypt` and `decrypt_all`.

Use ```cargo run --example endtoend``` to check correctness of the implementation.

**WARNING:** This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

## Dependencies
* [arkworks](http://arkworks.rs) project for finite field and elliptic curve arithmetic.
* [merlin](https://github.com/dalek-cryptography/merlin) library for implementing the Fiat-Shamir transform.

## Overview
* [`src/dealer`](src/dealer.rs): Contains an implementation of the `setup` and `epoch_setup` methods executed by a trusted dealer for the batched threshold encryption scheme.
* [`src/encryption`](src/encryption.rs): Contains an implementation of the `encrypt` method for the batched threshold encryption scheme.
* [`src/decryption`](src/decryption.rs): Contains an implementation of:
  * `partial_decrypt` - computed the message sent by each member of the committee to decrypt a batch of ciphertexts.
  * `decrypt_all` - gathers partial decryptions and recovers all messages from the batch of ciphertexts. This uses [FK22](https://github.com/khovratovich/Kate/blob/master/Kate_amortized.pdf) to compute all KZG opening proofs in quasi-linear time.

## License
This library is released under the MIT License.
