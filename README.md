# Matrix-based synchronous encryption method

A simple and reliable data encryption algorithm written in C++.

---

## ✨ Features

- Matrix-based.
- Synchronous encryption method.
- Simple and clean C++ API.
- Extended codebase, Gaussian method and matrix multiply.
- Easy to encrypt and decrypt data.
- Reliability algorithm.

---

## 📦 Getting Started

### 1. Clone the repository
```bash
git clone https://github.com/neitoroxin/matrixencrypt.git
cd matrixencrypt
```

---

## How use
input.txt - source message, string
output.txt - encrypted message
source.txt - decrypted message
encoder(/.exe) - code cipher text from the input file from input.txt to output.txt
decoder(/.exe) - program for decrypting messages from output.txt to source.txt

key.txt - private key for synchronous encryption, any vector
### hint: save space

---

## About algrotihm
The algorithm boils down to solving systems of linear algebraic equations. The key idea is that after solving a system, we cannot use its roots to determine the original coefficient matrix.
The presented implementation encrypts strings, but it can be generalized to any data arrays.

---

## More info (math and linear algebra):
let A - source message, B - private key. then solve equation Ax = B.
if we know roots, we cant restore source message from formula, because the matrices have different dimensions.
let M = B * x^T * (x^T * x)^{-1}. Then we can give as encrypted message A - M.

As for the algorithm's cryptographic strength, everything is excellent. The encryption method is synchronous and predicting a vector from a matrix A - M.
