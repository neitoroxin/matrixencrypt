#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <random>
#include <climits>
#include <filesystem>
#include <cstring>

class modular {
private:
    int value;
    static const int MOD = 257;

    int mod(int x) const {
        x %= MOD;
        return x < 0 ? x + MOD : x;
    }

public:
    modular(int x = 0) : value(mod(x)) {}

    operator int() const { return value; }

    modular operator+(int x) const { return modular(value + x); }
    modular operator-(int x) const { return modular(value - x); }
    modular operator*(int x) const { return modular((long long)value * x % MOD); }

    modular operator+(const modular& other) const { return *this + other.value; }
    modular operator-(const modular& other) const { return *this - other.value; }
    modular operator*(const modular& other) const { return *this * other.value; }

    modular& operator+=(int x) { value = mod(value + x); return *this; }
    modular& operator-=(int x) { value = mod(value - x); return *this; }
    modular& operator*=(int x) { value = mod((long long)value * x % MOD); return *this; }

    friend modular binpow(modular a, int k) {
        modular ans = 1;
        while (k) {
            if (k & 1) {
                ans *= a;
            }
            a *= a;
            k >>= 1;
        }
        return ans;
    }

    modular operator/(int x) const { return modular(value * binpow(modular(x), MOD - 2)); }
    modular operator/(const modular& other) const { return *this / other.value; }
    modular& operator/=(int x) { value = mod(modular(value * binpow(modular(x), MOD - 2))); return *this; }


    friend std::ostream& operator<<(std::ostream& os, const modular& m) {
        return os << m.value;
    }

    friend std::istream& operator>>(std::istream& is, modular& m) {
        int temp;
        is >> temp;
        m.value = m.mod(temp);
        return is;
    }
};

template <typename T>
struct matrix {
    std::vector<std::vector<T>> data;

    matrix(long long n, long long m);
    long long rows() const { return data.size(); }
    long long cols() const { return data[0].size(); }
    
    matrix(const matrix<T>& other) : data(other.data) {}
    matrix<T>& operator=(const matrix<T>& other) {
        if (this != &other) {
            data = other.data;
        }
        return *this;
    }
};

template <typename T>
matrix<T>::matrix(long long n, long long m) {
    data = std::vector<std::vector<T>>(n, std::vector<T>(m, 0));
}

template <typename T>
matrix<T> operator * (const matrix<T>& left, const matrix<T>& right) {
    matrix<T> mid(left.rows(), right.cols());
    for (long long i = 0; i < mid.rows(); i++) {
        for (long long j = 0; j < mid.cols(); j++) {
            for (long long r = 0; r < left.cols(); r++) {
                mid.data[i][j] += left.data[i][r] * right.data[r][j];
            }
        }
    }
    return mid;
}

template <typename T>
matrix<T> operator + (const matrix<T>& left, const matrix<T>& right) {
    matrix<T> mid(left.rows(), right.cols());
    for (long long i = 0; i < mid.rows(); i++) {
        for (long long j = 0; j < mid.cols(); j++) {
            mid.data[i][j] = left.data[i][j] + right.data[i][j];
        }
    }
    return mid;
}

template <typename T>
matrix<T> operator - (const matrix<T>& left, const matrix<T>& right) {
    matrix<T> mid(left.rows(), right.cols());
    for (long long i = 0; i < mid.rows(); i++) {
        for (long long j = 0; j < mid.cols(); j++) {
            mid.data[i][j] = left.data[i][j] - right.data[i][j];
        }
    }
    return mid;
}

template <typename T>
matrix<T> operator * (const T left, const matrix<T>& right) {
    matrix<T> mid(right.rows(), right.cols());
    for (long long i = 0; i < mid.rows(); i++) {
        for (long long j = 0; j < mid.cols(); j++) {
            mid.data[i][j] = right.data[i][j] * left;
        }
    }
    return mid;
}

template <typename T>
matrix<T> matrix_inverse(const matrix<T>& A) {
    const long long n = A.cols();
    matrix<T> A0(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A0.data[i][j] = modular(A.data[i][j]);
        }
    }

    matrix<T> result(n, n);
    for (int i = 0; i < n; i++) {
        result.data[i][i] = 1;
    }

    for (long long i = 0; i < n; i++) {
        long long best = i;
        for (long long j = i + 1; j < n; j++) {
            if (std::abs(A0.data[j][i]) > std::abs(A0.data[best][i])) {
                best = j;
            }
        }

        if (A0.data[best][i] == 0) {
            throw std::runtime_error("Matrix is singular, cannot invert");
        }

        std::swap(A0.data[i], A0.data[best]);
        std::swap(result.data[i], result.data[best]);

        T div = A0.data[i][i];
        for (long long j = 0; j < n; j++) {
            A0.data[i][j] /= div;
            result.data[i][j] /= div;
        }

        for (long long j = 0; j < n; j++) {
            if (j != i) {
                T fact = A0.data[j][i];
                for (long long k = 0; k < n; k++) {
                    A0.data[j][k] -= fact * A0.data[i][k];
                    result.data[j][k] -= fact * result.data[i][k];
                }
            }
        }
    }

    return result;
}

template <typename T>
matrix<T> transp(const matrix<T>& A) {
    long long rows = A.rows();
    long long cols = A.cols();
    matrix<T> result(cols, rows);

    for (long long i = 0; i < rows; i++) {
        for (long long j = 0; j < cols; j++) {
            result.data[j][i] = A.data[i][j];
        }
    }

    return result;
}

template <typename T>
matrix<T> gauss(const matrix<T>& A, const matrix<T>& B) {
    long long n = A.rows();
    long long m = B.cols();
    matrix<T> A0 = A;
    matrix<T> B0 = B;

    for (long long i = 0; i < n; i++) {
        long long best = i;
        for (long long j = i + 1; j < n; j++) {
            if (std::abs(A0.data[j][i]) > std::abs(A0.data[best][i])) {
                best = j;
            }
        }

        if (A0.data[best][i] == 0) {
            throw std::runtime_error("Matrix is singular in Gaussian elimination");
        }

        if (best != i) {
            std::swap(A0.data[i], A0.data[best]);
            std::swap(B0.data[i], B0.data[best]);
        }

        for (long long j = i + 1; j < n; j++) {
            T factor = A0.data[j][i] / A0.data[i][i];
            for (long long k = i; k < n; k++) {
                A0.data[j][k] -= factor * A0.data[i][k];
            }
            for (long long col = 0; col < m; col++) {
                B0.data[j][col] -= factor * B0.data[i][col];
            }
        }
    }

    matrix<T> result(n, m);
    for (long long i = n - 1; i >= 0; i--) {
        for (long long j = 0; j < m; j++) {
            T sum = 0;
            for (long long col = i + 1; col < n; col++) {
                sum += A0.data[i][col] * result.data[col][j];
            }
            result.data[i][j] = (T(B0.data[i][j]) - T(sum)) / T(A0.data[i][i]);
        }
    }

    return result;
}

template <typename T>
std::istream & operator >> (std::istream &in, matrix<T>& a) {
    for (long long i = 0; i < a.rows(); i++) {
        for (long long j = 0; j < a.cols(); j++) {
            in >> a.data[i][j];
        }
    }
    return in;
}

template <typename T>
std::ostream & operator << (std::ostream & out, const matrix<T>& a) {
    for (long long i = 0; i < a.rows(); i++) {
        for (long long j = 0; j < a.cols(); j++) {
            out << a.data[i][j] << ' ';
        }
        out << std::endl;
    }
    return out;
}

template <typename T>
matrix<T> fromTtoModular(const matrix<T>& A) {
    matrix<T> A0(A.rows(), A.cols());
    for (int i = 0; i < A0.rows(); i++) {
        for (int j = 0; j < A0.cols(); j++) {
            A0.data[i][j] = T(A.data[i][j]);
        }
    }
    return A0;
}

class Solution {
private:
    std::mt19937 mt;
    
public:
    Solution() : mt(std::random_device{}()) {}
    
    matrix<modular> maybeA(const matrix<modular>& B, const matrix<modular>& x) {
        const matrix<modular> xt = transp(x);
        matrix<modular> scalar1 = xt * x;
        
        if (scalar1.data[0][0] == 0) {
            throw std::runtime_error("Division by zero in maybeA calculation");
        }
        
        const modular scalar = 1 / scalar1.data[0][0];
        return scalar * (B * xt);
    }

    void encrypt(const std::string& input_file, const std::string& key_file, const std::string& output_file) {
        std::ifstream fin(input_file);
        std::ofstream fout(output_file);
        
        if (!fin.is_open()) {
            throw std::runtime_error("Cannot open input file: " + input_file);
        }
        if (!fout.is_open()) {
            throw std::runtime_error("Cannot open output file: " + output_file);
        }

        std::string s;
        std::getline(fin, s);
        long long n = s.size();
        long long matn = std::ceil(std::sqrt(n));

        std::uniform_int_distribution<char> dist(CHAR_MIN, CHAR_MAX);

        matrix<modular> A(matn, matn), B(matn, 1);
        
        // Filling matrix A
        for (int i = 0; i < n; i++) {
            A.data[i / matn][i % matn] = s[i] + 1;
        }

        // Filling random numbers
        for (int i = n; i < matn * matn; i++) {
            A.data[i / matn][i % matn] = dist(mt);
        }

        // Reading key
        std::ifstream key_stream(key_file);
        if (!key_stream.is_open()) {
            throw std::runtime_error("Cannot open key file: " + key_file);
        }
        
        for (int i = 0; i < matn; i++) {
            if (!(key_stream >> B.data[i][0])) {
                throw std::runtime_error("Invalid key file format");
            }
        }
        key_stream.close();

        // Encrypting
        const matrix<modular> x = gauss(A, B);
        matrix<modular> P = maybeA(B, x) - A;
        
        // Writing results
        fout << n << std::endl << transp(x) << std::endl << P << std::endl;
        
        fin.close();
        fout.close();
    }

    void decrypt(const std::string& encrypted_file, const std::string& key_file, const std::string& output_file) {
        std::ifstream fin(encrypted_file);
        std::ofstream fout(output_file);
        
        if (!fin.is_open()) {
            throw std::runtime_error("Cannot open encrypted file: " + encrypted_file);
        }
        if (!fout.is_open()) {
            throw std::runtime_error("Cannot open output file: " + output_file);
        }

        long long n;
        fin >> n;
        long long matn = std::ceil(std::sqrt(n));

        matrix<modular> x(matn, 1);
        for (int i = 0; i < matn; i++) {
            fin >> x.data[i][0];
        }

        // Reading key
        matrix<modular> B(matn, 1);
        std::ifstream key_stream(key_file);
        if (!key_stream.is_open()) {
            throw std::runtime_error("Cannot open key file: " + key_file);
        }
        
        for (int i = 0; i < matn; i++) {
            if (!(key_stream >> B.data[i][0])) {
                throw std::runtime_error("Invalid key file format");
            }
        }
        key_stream.close();

        const matrix<modular> maybe = maybeA(B, x);
        
        // Reading and encrypting
        for (int i = 0; i < matn; i++) {
            for (int j = 0; j < matn; j++) {
                modular reserv;
                if (!(fin >> reserv)) {
                    throw std::runtime_error("Invalid encrypted file format");
                }
                if (i * matn + j < n) {
                    char decrypted_char = static_cast<char>(maybe.data[i][j] - reserv - 1);
                    fout << decrypted_char;
                }
            }
        }
        
        fin.close();
        fout.close();
    }
};

void print_usage() {
    std::cout << "Usage:\n"
              << "  matrixencrypt -encrypt <key_file> -o <output_file> <input_file>\n"
              << "  matrixencrypt -decrypt <key_file> -o <output_file> <encrypted_file>\n"
              << "Examples:\n"
              << "  matrixencrypt -encrypt key.txt -o encrypted.txt message.txt\n"
              << "  matrixencrypt -decrypt key.txt -o decrypted.txt encrypted.txt\n";
}

int main(int argc, char* argv[]) {

    if (argc < 6) {
        print_usage();
        return 1;
    }

    try {
        Solution crypto;
        std::string mode = argv[1];
        std::string key_file = argv[2];
        std::string output_flag = argv[3];
        std::string output_file = argv[4];
        std::string input_file = argv[5];

        if (output_flag != "-o") {
            std::cerr << "expected '-o' flag before output file\n";
            print_usage();
            return 1;
        }

        if (mode == "-encrypt") {
            crypto.encrypt(input_file, key_file, output_file);
            std::cout << "File encrypted successfully: " << output_file << std::endl;
        } else if (mode == "-decrypt") {
            crypto.decrypt(input_file, key_file, output_file);
            std::cout << "File decrypted successfully: " << output_file << std::endl;
        } else {
            std::cerr << "unknown mode " << mode << "'\n";
            print_usage();
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
