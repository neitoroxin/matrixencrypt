#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <random>
#include <climits>

template <typename T>
struct matrix {
    std::vector<std::vector<T>> data;

    matrix(long long n, long long m);
    long long rows() const {
		return data.size();
	}
    long long cols() const {
		return data[0].size();
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
matrix<double> pow(const matrix<T>& A, const T pow) {
    const long long n = A.cols();
    matrix<double> A0(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A0.data[i][j] = double(A.data[i][j]);
        }
    }

    if (pow != -1) return A0;
    matrix<double> result(n, n);

    for (int i = 0; i < n; i++) {
        result.data[i][i] = 1;
    }

    for (long long i = 0; i < n; i++) {
        long long best = i;
        for (long long j = i + 1; j < n; j++) {
            if (abs(A0.data[j][i]) > abs(A.data[best][i])) {
                best = j;
            }
        }

        std::swap(A0.data[i], A0.data[best]);
        std::swap(result.data[i], result.data[best]);

        double div = A0.data[i][i];
        for (long long j = 0; j < n; j++) {
            A0.data[i][j] /= div;
            result.data[i][j] /= div;
        }

        for (long long j = 0; j < n; j++) {
            if (j != i) {
                double fact = A0.data[j][i];
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
matrix<double> gauss(const matrix<T>& A, const matrix<T>& B) {
    long long n = A.rows();
    long long m = B.cols();
    matrix<T> A0 = A;
    matrix<T> B0 = B;

    for (long long i = 0; i < n; i++) {
        long long best = i;
        for (long long j = i + 1; j < n; j++) {
            if (abs(A0.data[j][i]) > abs(A0.data[best][i])) {
                best = j;
            }
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

    matrix<double> result(n, m);
    for (long long i = n - 1; i >= 0; i--) {
        for (long long j = 0; j < m; j++) {
            T sum = 0;
            for (long long col = i + 1; col < n; col++) {
                sum += A0.data[i][col] * result.data[col][j];
            }
            result.data[i][j] =  (double(B0.data[i][j]) - double(sum)) / double(A0.data[i][i]);
        }
    }

    return result;
}

template <typename T>
std::ostream & operator << (std::ostream & out, const matrix<T>& a) { // Вывод матрицы
    for (long long i = 0; i < a.rows(); i++) {
        for (long long j = 0; j < a.cols(); j++) {
            out << a.data[i][j] << ' ';
        }
        out << '\n';
    }
    return out;
}

template <typename T>
matrix<double> fromTtoDouble(const matrix<T>& A) {
    matrix<double> A0(A.rows(), A.cols());
    for (int i = 0; i < A0.rows(); i++) {
        for (int j = 0; j < A0.cols(); j++) {
            A0.data[i][j] = double(A.data[i][j]);
        }
    }

    return A0;
}


std::string s;
long long n, matn;
std::mt19937 mt(time(nullptr));
std::uniform_int_distribution<char> dist(CHAR_MIN, CHAR_MAX);

class Solution {
public:
    Solution() {}

    template <typename T>
    matrix<T> maybeA(const matrix<T>& B, const matrix<T>& x) {
        const matrix<T> xt = transp(x);
        matrix<double> scalar1 = xt * x;
        const double scalar = 1 / scalar1.data[0][0];

        return scalar * (B * xt);
    }

    void encoder() { // Encoder with string
        std::ifstream fin("input.txt"); // File with source message
        std::ofstream fout("output.txt"); // File with encrypted message
        getline(fin, s);
        n = s.size();
        matn = ceil(sqrt(n));

        matrix<long long> A(matn, matn), B(matn, 1);
        for (int i = 0; i < n; i++) {
            A.data[i / matn][i % matn] = s[i] + 1;
        }

        for (int i = n; i < matn * matn; i++) {
            A.data[i / matn][i % matn] = dist(mt);
        }

        std::cerr << A << std::endl;

        std::ifstream is("key.txt"); // Reading key
        for (int i = 0; i < matn; i++) {
            is >> B.data[i][0];
        }
        is.close();

        const matrix<double> x{gauss(A, B)};
        fout << n << std::endl << transp(x) << std::endl << maybeA(fromTtoDouble(B), x) - fromTtoDouble(A) << std::endl;
        fin.close();
        fout.close();
    }

    void decoder() {
        std::ifstream fin("output.txt"); // File with encrypted message
        std::ofstream fout("source.txt"); // File with decrypted message
        fin >> n;
        matn = ceil(sqrt(n));

        matrix<double> x(matn, 1);
        for (int i = 0; i < matn; i++) {
            fin >> x.data[i][0];
        }
        matrix<double> xt = transp(x);
        matrix<long long> B(matn, 1);

        std::ifstream is("key.txt"); // Reading key
        for (int i = 0; i < matn; i++) {
            is >> B.data[i][0];
        }
        is.close();

        const matrix<double> maybe{maybeA(fromTtoDouble(B), x)};
        for (int i = 0; i < matn; i++) {
            for (int j = 0; j < matn; j++) {
                double reserv;
                fin >> reserv;
                if (i * matn + j < n) {
                    fout << char(round(maybe.data[i][j] - reserv - 1));
                } else {
                    fin.close();
                    fout.close();
                    return;
                }
            }
        }
        fin.close();
        fout.close();
    }
};

#include <chrono>
int main() {
    std::ios::sync_with_stdio(0);
    std::cin.tie(0);
    std::cout.tie(0);

    auto start = std::chrono::high_resolution_clock::now();
    Solution s = Solution();

    std::cout << std::endl;

    int flag = 0;
    if (flag == 1) {
        s.encoder();
    } else if (flag == 0) {
        s.decoder();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << "\n\n" << std::chrono::duration<double, std::milli>(end - start).count() << " ms\n";
    return 0;
}