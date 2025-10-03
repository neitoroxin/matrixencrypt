# MatrixEncrypt

MatrixEncrypt is a C++ project that demonstrates **synchronous encryption using linear algebra**.  
The algorithm is based on solving matrix equations of the form:


where:
- `A` is the encryption key (matrix),
- `x` is the plaintext vector (input message),
- `B` is the ciphertext vector (output message).

The project provides two programs:  
- `encoder` — encrypts a message using a key matrix  
- `decoder` — decrypts a ciphertext back into the original message  

---

## ✨ Features

- Encrypts text data using matrix operations  
- Decrypts ciphertext back to original plaintext  
- Works with text files as input/output  
- Simple implementation in C++ for educational purposes  

---

## 📦 Getting Started

### 1. Clone the repository
```bash
git clone https://github.com/neitoroxin/matrixencrypt.git
cd matrixencrypt
```

## How use

### All programs must be located in one directory.

- Make sure all programs are in the same folder.
- When you write you message in 'input.txt', run 'encoder'. Encrypted message will appear in 'output.txt'
- For decrypt encrypted message use 'decoder', when encrypted message will be in 'output.txt'. The decrypted message will appear in 'source.txt'
