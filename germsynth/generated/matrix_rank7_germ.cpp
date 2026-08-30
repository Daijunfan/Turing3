// Generated from matrix_rank7_germ.json; do not edit coefficient arrays by hand.
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

using Scalar = long long;
using Matrix = std::vector<Scalar>;
constexpr int RANK = 7;
constexpr std::array<std::array<int, 4>, RANK> U = {{
    {1, 0, 0, 0},
    {1, 1, 0, 0},
    {1, 0, -1, 0},
    {0, 0, 0, 1},
    {1, 0, 0, 1},
    {0, 1, 0, -1},
    {0, 0, 1, 1}
}};
constexpr std::array<std::array<int, 4>, RANK> V = {{
    {0, 1, 0, -1},
    {0, 0, 0, 1},
    {1, 1, 0, 0},
    {1, 0, -1, 0},
    {1, 0, 0, 1},
    {0, 0, 1, 1},
    {1, 0, 0, 0}
}};
constexpr std::array<std::array<int, 4>, RANK> W = {{
    {0, 1, 0, 1},
    {-1, 1, 0, 0},
    {0, 0, 0, -1},
    {-1, 0, -1, 0},
    {1, 0, 0, 1},
    {1, 0, 0, 0},
    {0, 0, 1, -1}
}};

struct Counter { std::uint64_t scalar_multiplications = 0; };

static Matrix zero_matrix(int n) { return Matrix(static_cast<std::size_t>(n) * n, 0); }
static Scalar& at(Matrix& a, int n, int row, int col) { return a[static_cast<std::size_t>(row) * n + col]; }
static Scalar at(const Matrix& a, int n, int row, int col) { return a[static_cast<std::size_t>(row) * n + col]; }

static std::array<Matrix, 4> split(const Matrix& a, int n) {
    const int h = n / 2;
    std::array<Matrix, 4> blocks = {zero_matrix(h), zero_matrix(h), zero_matrix(h), zero_matrix(h)};
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < h; ++j) {
            at(blocks[0], h, i, j) = at(a, n, i, j);
            at(blocks[1], h, i, j) = at(a, n, i, j + h);
            at(blocks[2], h, i, j) = at(a, n, i + h, j);
            at(blocks[3], h, i, j) = at(a, n, i + h, j + h);
        }
    }
    return blocks;
}

static Matrix join(const std::array<Matrix, 4>& blocks, int h) {
    Matrix result = zero_matrix(2 * h);
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < h; ++j) {
            at(result, 2 * h, i, j) = at(blocks[0], h, i, j);
            at(result, 2 * h, i, j + h) = at(blocks[1], h, i, j);
            at(result, 2 * h, i + h, j) = at(blocks[2], h, i, j);
            at(result, 2 * h, i + h, j + h) = at(blocks[3], h, i, j);
        }
    }
    return result;
}

static Matrix linear_form(const std::array<int, 4>& coefficients,
                          const std::array<Matrix, 4>& blocks, int n) {
    Matrix result = zero_matrix(n);
    for (int block = 0; block < 4; ++block) {
        if (coefficients[block] == 0) continue;
        for (std::size_t index = 0; index < result.size(); ++index) {
            result[index] += static_cast<Scalar>(coefficients[block]) * blocks[block][index];
        }
    }
    return result;
}

static Matrix germ_multiply(const Matrix& a, const Matrix& b, int n, Counter& counter) {
    if (n == 1) {
        ++counter.scalar_multiplications;
        return Matrix{a[0] * b[0]};
    }
    const int h = n / 2;
    const auto a_blocks = split(a, n);
    const auto b_blocks = split(b, n);
    std::array<Matrix, RANK> products;
    for (int term = 0; term < RANK; ++term) {
        Matrix left = linear_form(U[term], a_blocks, h);
        Matrix right = linear_form(V[term], b_blocks, h);
        products[term] = germ_multiply(left, right, h, counter);
    }
    std::array<Matrix, 4> output = {zero_matrix(h), zero_matrix(h), zero_matrix(h), zero_matrix(h)};
    for (int term = 0; term < RANK; ++term) {
        for (int block = 0; block < 4; ++block) {
            if (W[term][block] == 0) continue;
            for (std::size_t index = 0; index < output[block].size(); ++index) {
                output[block][index] += static_cast<Scalar>(W[term][block]) * products[term][index];
            }
        }
    }
    return join(output, h);
}

static Matrix schoolbook(const Matrix& a, const Matrix& b, int n) {
    Matrix result = zero_matrix(n);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < n; ++j)
                at(result, n, i, j) += at(a, n, i, k) * at(b, n, k, j);
    return result;
}

static bool verify_tensor_identity() {
    std::array<int, 64> reconstructed{};
    std::array<int, 64> target{};
    for (int term = 0; term < RANK; ++term)
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                    reconstructed[(i * 4 + j) * 4 + k] += U[term][i] * V[term][j] * W[term][k];
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k) {
                const int ai = i * 2 + k;
                const int bj = k * 2 + j;
                const int ci = i * 2 + j;
                ++target[(ai * 4 + bj) * 4 + ci];
            }
    return reconstructed == target;
}

static std::uint64_t integer_power(std::uint64_t base, int exponent) {
    std::uint64_t value = 1;
    for (int i = 0; i < exponent; ++i) value *= base;
    return value;
}

int main() {
    if (!verify_tensor_identity()) {
        std::cerr << "tensor identity: FAIL\n";
        return 1;
    }
    std::mt19937_64 rng(20260830ULL);
    std::uniform_int_distribution<int> distribution(-2, 2);
    for (int n : {1, 2, 4, 8, 16, 32, 64, 128}) {
        Matrix a = zero_matrix(n), b = zero_matrix(n);
        for (auto& value : a) value = distribution(rng);
        for (auto& value : b) value = distribution(rng);
        Counter counter;
        Matrix actual = germ_multiply(a, b, n, counter);
        Matrix expected = schoolbook(a, b, n);
        const int depth = [] (int x) { int d = 0; while (x > 1) { x /= 2; ++d; } return d; }(n);
        const std::uint64_t predicted = integer_power(RANK, depth);
        if (actual != expected || counter.scalar_multiplications != predicted) {
            std::cerr << "recursive verification failed at n=" << n << "\n";
            return 2;
        }
        std::cout << "n=" << n << " exact=true scalar_multiplications="
                  << counter.scalar_multiplications << "\n";
    }
    std::cout << "CXX_INDEPENDENT_VERIFICATION=PASS\n";
    return 0;
}
