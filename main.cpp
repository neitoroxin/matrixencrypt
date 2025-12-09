#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <random>
#include <climits>
#include <ctime>
#include <filesystem>
#include <cstring>
#include <chrono>

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

    modular operator/(int x) const { 
        if (x == 0) throw std::runtime_error("Division by zero");
        return *this * modular(x).inverse();
    }
    
    modular operator/(const modular& other) const { 
        return *this * other.inverse();
    }
    
    modular& operator/=(int x) { 
        if (x == 0) throw std::runtime_error("Division by zero");
        *this *= modular(x).inverse();
        return *this;
    }

    modular inverse() const {
        return pow(MOD - 2);
    }
    
    modular pow(int k) const {
        modular result = 1;
        modular base = *this;
        while (k > 0) {
            if (k & 1) result *= base;
            base *= base;
            k >>= 1;
        }
        return result;
    }
    
    bool operator==(const modular& other) const { return value == other.value; }
    bool operator!=(const modular& other) const { return value != other.value; }
    bool operator<(const modular& other) const { return value < other.value; }
    bool operator>(const modular& other) const { return value > other.value; }

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

    matrix(long long n = 0, long long m = 0) : data(n, std::vector<T>(m, 0)) {}
    
    long long rows() const { return data.size(); }
    long long cols() const { return data.empty() ? 0 : data[0].size(); }
    
    matrix(const matrix<T>& other) = default;
    matrix<T>& operator=(const matrix<T>& other) = default;
    
    std::vector<T>& operator[](size_t i) { return data[i]; }
    const std::vector<T>& operator[](size_t i) const { return data[i]; }
    
    void resize(long long n, long long m) {
        data.resize(n, std::vector<T>(m, 0));
    }
};

template <typename T>
matrix<T> operator*(const matrix<T>& left, const matrix<T>& right) {
    if (left.cols() != right.rows()) {
        throw std::runtime_error("Matrix dimensions mismatch for multiplication");
    }
    
    matrix<T> result(left.rows(), right.cols());
    for (long long i = 0; i < result.rows(); i++) {
        for (long long j = 0; j < result.cols(); j++) {
            T sum = 0;
            for (long long k = 0; k < left.cols(); k++) {
                sum += left[i][k] * right[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

template <typename T>
matrix<T> operator+(const matrix<T>& left, const matrix<T>& right) {
    if (left.rows() != right.rows() || left.cols() != right.cols()) {
        throw std::runtime_error("Matrix dimensions mismatch for addition");
    }
    
    matrix<T> result(left.rows(), left.cols());
    for (long long i = 0; i < result.rows(); i++) {
        for (long long j = 0; j < result.cols(); j++) {
            result[i][j] = left[i][j] + right[i][j];
        }
    }
    return result;
}

template <typename T>
matrix<T> operator-(const matrix<T>& left, const matrix<T>& right) {
    if (left.rows() != right.rows() || left.cols() != right.cols()) {
        throw std::runtime_error("Matrix dimensions mismatch for subtraction");
    }
    
    matrix<T> result(left.rows(), left.cols());
    for (long long i = 0; i < result.rows(); i++) {
        for (long long j = 0; j < result.cols(); j++) {
            result[i][j] = left[i][j] - right[i][j];
        }
    }
    return result;
}

template <typename T>
matrix<T> operator*(const T& scalar, const matrix<T>& mat) {
    matrix<T> result(mat.rows(), mat.cols());
    for (long long i = 0; i < result.rows(); i++) {
        for (long long j = 0; j < result.cols(); j++) {
            result[i][j] = scalar * mat[i][j];
        }
    }
    return result;
}

template <typename T>
matrix<T> matrix_inverse(const matrix<T>& A) {
    long long n = A.rows();
    if (n != A.cols()) {
        throw std::runtime_error("Matrix must be square for inversion");
    }
    
    matrix<T> aug(n, 2 * n);
    for (long long i = 0; i < n; i++) {
        for (long long j = 0; j < n; j++) {
            aug[i][j] = A[i][j];
        }
        aug[i][n + i] = 1;
    }
    
    for (long long i = 0; i < n; i++) {
        if (aug[i][i] == modular(0)) {
            long long swap_row = -1;
            for (long long j = i + 1; j < n; j++) {
                if (aug[j][i] != modular(0)) {
                    swap_row = j;
                    break;
                }
            }
            if (swap_row == -1) {
                throw std::runtime_error("Matrix is singular, cannot invert");
            }
            std::swap(aug[i], aug[swap_row]);
        }
        
        T pivot = aug[i][i];
        for (long long j = 0; j < 2 * n; j++) {
            aug[i][j] /= pivot;
        }
        
        for (long long j = 0; j < n; j++) {
            if (j != i) {
                T factor = aug[j][i];
                for (long long k = 0; k < 2 * n; k++) {
                    aug[j][k] -= factor * aug[i][k];
                }
            }
        }
    }
    
    matrix<T> inv(n, n);
    for (long long i = 0; i < n; i++) {
        for (long long j = 0; j < n; j++) {
            inv[i][j] = aug[i][n + j];
        }
    }
    
    return inv;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const matrix<T>& a) {
    for (long long i = 0; i < a.rows(); i++) {
        for (long long j = 0; j < a.cols(); j++) {
            out << a[i][j];
            if (j < a.cols() - 1) out << ' ';
        }
        if (i < a.rows() - 1) out << '\n';
    }
    return out;
}

template <typename T>
std::istream& operator>>(std::istream& in, matrix<T>& a) {
    for (long long i = 0; i < a.rows(); i++) {
        for (long long j = 0; j < a.cols(); j++) {
            if (!(in >> a[i][j])) {
                throw std::runtime_error("Failed to read matrix element");
            }
        }
    }
    return in;
}

class Solution {
private:
    std::mt19937 rng;
    const int MOD = 257;
    
    int random_int() {
        std::uniform_int_distribution<int> dist(1, MOD - 2);
        return dist(rng);
    }
    
    matrix<modular> generate_random_matrix(int n) {
        matrix<modular> result(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = modular(random_int());
            }
        }
        return result;
    }
    
    bool commute(const matrix<modular>& A, const matrix<modular>& B) {
        matrix<modular> AB = A * B;
        matrix<modular> BA = B * A;
        
        for (int i = 0; i < A.rows(); i++) {
            for (int j = 0; j < A.cols(); j++) {
                if (AB[i][j] != BA[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }
    
    matrix<modular> matrix_power(const matrix<modular>& M, int power) {
        int n = M.rows();
        matrix<modular> result(n, n);
        matrix<modular> base = M;
        
        // init once matrix
        for (int i = 0; i < n; i++) {
            result[i][i] = modular(1);
        }
        
        if (power < 0) {
            base = matrix_inverse(M);
            power = -power;
        }
        
        while (power > 0) {
            if (power & 1) {
                result = result * base;
            }
            base = base * base;
            power >>= 1;
        }
        
        return result;
    }
    
public:
    Solution() : rng(std::chrono::steady_clock::now().time_since_epoch().count()) {}

    void generate_key(int n, const std::string& key_file) {
        std::ofstream fout(key_file);
        if (!fout.is_open()) {
            throw std::runtime_error("Cannot open output file: " + key_file);
        }

        matrix<modular> G, H;
        
        while (true) {
            matrix<modular> P = generate_random_matrix(n);
            
            // generate diagonal matrix with power 3 and 5
            matrix<modular> D1(n, n), D2(n, n);
            for (int i = 0; i < n; i++) {
                D1[i][i] = modular(3).pow((i + 1) * 7 % (MOD - 1));
                D2[i][i] = modular(5).pow((i + 1) * 11 % (MOD - 1));
            }
            
            // calcucation comutate matrix
            matrix<modular> P_inv = matrix_inverse(P);
            G = P * D1 * P_inv;
            H = P * D2 * P_inv;
            
            if (commute(G, H)) break;
        }
        
        int a = random_int();
        int b = random_int();
        
        // calcucation open keys with binary powing
        matrix<modular> P_A = matrix_power(G, a);
        matrix<modular> P_B = matrix_power(H, b);
        
        fout << n << std::endl;
        fout << G << std::endl;
        fout << H << std::endl;
        fout << P_A << std::endl;
        fout << P_B << std::endl;
        fout << a << " " << b << std::endl;
        
        fout.close();
    }

    void encrypt(const std::string& input_file, const std::string& key_file, 
                const std::string& output_file) {
        std::ifstream fin(input_file);
        std::ifstream kin(key_file);
        std::ofstream fout(output_file);
        
        if (!fin.is_open()) throw std::runtime_error("Cannot open input file: " + input_file);
        if (!kin.is_open()) throw std::runtime_error("Cannot open key file: " + key_file);
        if (!fout.is_open()) throw std::runtime_error("Cannot open output file: " + output_file);

        int n;
        kin >> n;
        
        matrix<modular> G(n, n), H(n, n), P_A(n, n), P_B(n, n);
        kin >> G >> H >> P_A >> P_B;
        
        int a, b;
        kin >> a >> b;
        
        int m;
        fin >> m;
        if (m != n) throw std::runtime_error("Matrix size mismatch");
        
        matrix<modular> M(n, n);
        fin >> M;
        
        int r = random_int();
        int s = random_int();
        
        // encrypting with binary powing
        matrix<modular> C1 = matrix_power(G, r);
        matrix<modular> C2 = matrix_power(H, s);
        matrix<modular> P_A_pow_r = matrix_power(P_A, r);
        matrix<modular> P_B_pow_s = matrix_power(P_B, s);
        matrix<modular> C3 = P_A_pow_r * M * P_B_pow_s;
        
        fout << n << std::endl;
        fout << C1 << std::endl;
        fout << C2 << std::endl;
        fout << C3 << std::endl;
        fout << r << " " << s << std::endl;
        
        fin.close();
        kin.close();
        fout.close();
    }

    void decrypt(const std::string& input_file, const std::string& key_file,
                const std::string& output_file) {
        std::ifstream fin(input_file);
        std::ifstream kin(key_file);
        std::ofstream fout(output_file);
        
        if (!fin.is_open()) throw std::runtime_error("Cannot open input file: " + input_file);
        if (!kin.is_open()) throw std::runtime_error("Cannot open key file: " + key_file);
        if (!fout.is_open()) throw std::runtime_error("Cannot open output file: " + output_file);

        int n;
        kin >> n;
        
        matrix<modular> G(n, n), H(n, n), P_A(n, n), P_B(n, n);
        kin >> G >> H >> P_A >> P_B;
        
        int a, b;
        kin >> a >> b;
        
        fin >> n;
        
        matrix<modular> C1(n, n), C2(n, n), C3(n, n);
        fin >> C1 >> C2 >> C3;
        
        int r, s;
        fin >> r >> s;
        
        // decrypting with using binary powing
        matrix<modular> C1_inv_a = matrix_power(C1, -a);
        matrix<modular> C2_inv_b = matrix_power(C2, -b);
        
        matrix<modular> M = C1_inv_a * C3 * C2_inv_b;
        
        fout << n << std::endl << M << std::endl;
        
        fin.close();
        kin.close();
        fout.close();
    }

    void homomorphic_multiply(const std::string& input1_file, const std::string& input2_file,
                             const std::string& output_file) {
        std::ifstream fin1(input1_file);
        std::ifstream fin2(input2_file);
        std::ofstream fout(output_file);
        
        if (!fin1.is_open()) throw std::runtime_error("Cannot open input file 1: " + input1_file);
        if (!fin2.is_open()) throw std::runtime_error("Cannot open input file 2: " + input2_file);
        if (!fout.is_open()) throw std::runtime_error("Cannot open output file: " + output_file);

        int n;
        fin1 >> n;
        
        matrix<modular> C1_1(n, n), C2_1(n, n), C3_1(n, n);
        matrix<modular> C1_2(n, n), C2_2(n, n), C3_2(n, n);
        
        fin1 >> C1_1 >> C2_1 >> C3_1;
        fin2 >> n >> C1_2 >> C2_2 >> C3_2;
        
        // homo multiply:
        matrix<modular> C1_new = C1_1 * C1_2;
        matrix<modular> C2_new = C2_1 * C2_2;
        matrix<modular> C3_new = C3_1 * C3_2;
        
        fout << n << std::endl;
        fout << C1_new << std::endl;
        fout << C2_new << std::endl;
        fout << C3_new << std::endl;
        
        fin1.close();
        fin2.close();
        fout.close();
    }

    void print_usage() {
        std::cout << "Usage:\n"
                  << "  ./matrixcrypto -genkey <n> -o <key_file>\n"
                  << "  ./matrixcrypto -encrypt <key_file> -o <output_file> <input_file>\n"
                  << "  ./matrixcrypto -decrypt <key_file> -o <output_file> <input_file>\n"
                  << "  ./matrixcrypto -hmult <ciphertext1> <ciphertext2> -o <output_file>\n";
    }

    int main(int argc, char* argv[]) {
        if (argc < 2) {
            print_usage();
            return 1;
        }

        try {
            std::string mode = argv[1];
            
            if (mode == "-genkey" && argc == 5) {
                int n = std::stoi(argv[2]);
                std::string output_flag = argv[3];
                std::string key_file = argv[4];
                
                if (output_flag != "-o") {
                    std::cerr << "Error: expected '-o' flag before output file\n";
                    print_usage();
                    return 1;
                }
                
                generate_key(n, key_file);
                
            } else if (mode == "-encrypt" && argc == 6) {
                std::string key_file = argv[2];
                std::string output_flag = argv[3];
                std::string output_file = argv[4];
                std::string input_file = argv[5];
                
                if (output_flag != "-o") {
                    std::cerr << "Error: expected '-o' flag before output file\n";
                    print_usage();
                    return 1;
                }
                
                encrypt(input_file, key_file, output_file);
                
            } else if (mode == "-decrypt" && argc == 6) {
                std::string key_file = argv[2];
                std::string output_flag = argv[3];
                std::string output_file = argv[4];
                std::string input_file = argv[5];
                
                if (output_flag != "-o") {
                    std::cerr << "Error: expected '-o' flag before output file\n";
                    print_usage();
                    return 1;
                }
                
                decrypt(input_file, key_file, output_file);
                
            } else if (mode == "-hmult" && argc == 5) {
                std::string input1_file = argv[2];
                std::string input2_file = argv[3];
                std::string output_file = argv[4];
                
                homomorphic_multiply(input1_file, input2_file, output_file);
                
            } else {
                std::cerr << "Error: invalid arguments\n";
                print_usage();
                return 1;
            }
            
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
        
        return 0;
    }
};

int main(int argc, char* argv[]) {
    Solution app;
    return app.main(argc, argv);
}