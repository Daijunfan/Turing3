#!/usr/bin/env python3
"""Emit a standalone C++17 verifier for a polymorphic germ-pool certificate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def cpp_vec(values: list[int]) -> str:
    return "{" + ", ".join(str(value) for value in values) + "}"


def cpp_terms(terms: list[dict[str, list[int]]]) -> str:
    rows = []
    for term in terms:
        rows.append("Term{" + cpp_vec(term["u"]) + ", " + cpp_vec(term["v"]) + ", " + cpp_vec(term["w"]) + "}")
    return "std::array<Term, RANK>{{\n        " + ",\n        ".join(rows) + "\n    }}"


def emit(certificate: Path, output: Path) -> None:
    data = json.loads(certificate.read_text(encoding="utf-8"))
    if data.get("schema") != "germsynth-phenotype-pool-v1":
        raise ValueError("unsupported certificate schema")
    phenotypes = data["phenotypes"]
    pool_size = len(phenotypes)
    rank = data["rank"]
    selected = data["minimum_single_fault_cover"]["phenotype_indices"]
    phenotype_initializers = ",\n    ".join(cpp_terms(item["terms"]) for item in phenotypes)
    selected_init = cpp_vec(selected)

    source = f'''// Generated from {certificate.name}. Verification logic is standalone C++17.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

using Scalar = long long;
using Matrix = std::vector<Scalar>;
constexpr int POOL = {pool_size};
constexpr int RANK = {rank};
constexpr int COVER = {len(selected)};

struct Term {{
    std::array<int, 4> u;
    std::array<int, 4> v;
    std::array<int, 4> w;
}};

constexpr std::array<std::array<Term, RANK>, POOL> GERMS = {{{{
    {phenotype_initializers}
}}}};
constexpr std::array<int, COVER> SELECTED = {selected_init};

struct Counter {{ std::uint64_t scalar_multiplications = 0; }};

static std::uint64_t splitmix64(std::uint64_t x) {{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}}

static int support_mask(const std::array<int, 4>& coefficients) {{
    int mask = 0;
    for (int i = 0; i < 4; ++i) if (coefficients[i] != 0) mask |= (1 << i);
    return mask;
}}

// Resources are encoded as mode*16 + support mask: mode 0=left,1=right,2=output.
static std::set<int> resources_of(int phenotype) {{
    std::set<int> result;
    for (const Term& term : GERMS[phenotype]) {{
        result.insert(0 * 16 + support_mask(term.u));
        result.insert(1 * 16 + support_mask(term.v));
        result.insert(2 * 16 + support_mask(term.w));
    }}
    return result;
}}

static bool survives(int phenotype, const std::vector<int>& failed) {{
    const auto required = resources_of(phenotype);
    for (int resource : failed) if (required.count(resource)) return false;
    return true;
}}

static bool any_survives(const std::vector<int>& selected, const std::vector<int>& failed) {{
    for (int phenotype : selected) if (survives(phenotype, failed)) return true;
    return false;
}}

static std::vector<int> resource_universe() {{
    std::set<int> universe;
    for (int p = 0; p < POOL; ++p) {{
        const auto required = resources_of(p);
        universe.insert(required.begin(), required.end());
    }}
    return std::vector<int>(universe.begin(), universe.end());
}}

static bool verify_tensor_identity(int phenotype) {{
    std::array<int, 64> reconstructed{{}};
    std::array<int, 64> target{{}};
    for (const Term& term : GERMS[phenotype])
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                    reconstructed[(i * 4 + j) * 4 + k] += term.u[i] * term.v[j] * term.w[k];
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k) {{
                const int ai = i * 2 + k;
                const int bj = k * 2 + j;
                const int ci = i * 2 + j;
                ++target[(ai * 4 + bj) * 4 + ci];
            }}
    return reconstructed == target;
}}

static Matrix zero_matrix(int n) {{ return Matrix(static_cast<std::size_t>(n) * n, 0); }}
static Scalar& at(Matrix& a, int n, int row, int col) {{ return a[static_cast<std::size_t>(row) * n + col]; }}
static Scalar at(const Matrix& a, int n, int row, int col) {{ return a[static_cast<std::size_t>(row) * n + col]; }}

static std::array<Matrix, 4> split(const Matrix& a, int n) {{
    const int h = n / 2;
    std::array<Matrix, 4> blocks = {{zero_matrix(h), zero_matrix(h), zero_matrix(h), zero_matrix(h)}};
    for (int i = 0; i < h; ++i)
        for (int j = 0; j < h; ++j) {{
            at(blocks[0], h, i, j) = at(a, n, i, j);
            at(blocks[1], h, i, j) = at(a, n, i, j + h);
            at(blocks[2], h, i, j) = at(a, n, i + h, j);
            at(blocks[3], h, i, j) = at(a, n, i + h, j + h);
        }}
    return blocks;
}}

static Matrix join(const std::array<Matrix, 4>& blocks, int h) {{
    Matrix result = zero_matrix(2 * h);
    for (int i = 0; i < h; ++i)
        for (int j = 0; j < h; ++j) {{
            at(result, 2 * h, i, j) = at(blocks[0], h, i, j);
            at(result, 2 * h, i, j + h) = at(blocks[1], h, i, j);
            at(result, 2 * h, i + h, j) = at(blocks[2], h, i, j);
            at(result, 2 * h, i + h, j + h) = at(blocks[3], h, i, j);
        }}
    return result;
}}

static Matrix linear_form(const std::array<int, 4>& coefficients,
                          const std::array<Matrix, 4>& blocks, int n) {{
    Matrix result = zero_matrix(n);
    for (int block = 0; block < 4; ++block) {{
        if (coefficients[block] == 0) continue;
        for (std::size_t i = 0; i < result.size(); ++i)
            result[i] += static_cast<Scalar>(coefficients[block]) * blocks[block][i];
    }}
    return result;
}}

enum class Mode {{ Mixed, FaultAware }};

static int choose_phenotype(Mode mode, std::uint64_t node_id, const std::vector<int>& universe) {{
    if (mode == Mode::Mixed) return static_cast<int>(splitmix64(node_id) % POOL);
    const int failed = universe[splitmix64(node_id ^ 0xD1B54A32D192ED03ULL) % universe.size()];
    for (int phenotype : SELECTED) if (survives(phenotype, std::vector<int>{{failed}})) return phenotype;
    throw std::runtime_error("minimum cover failed to survive a single resource fault");
}}

static Matrix polymorphic_multiply(const Matrix& a, const Matrix& b, int n,
                                   Counter& counter, Mode mode,
                                   const std::vector<int>& universe,
                                   std::uint64_t node_id = 0x123456789abcdefULL) {{
    if (n == 1) {{
        ++counter.scalar_multiplications;
        return Matrix{{a[0] * b[0]}};
    }}
    const int phenotype = choose_phenotype(mode, node_id, universe);
    const int h = n / 2;
    const auto ab = split(a, n);
    const auto bb = split(b, n);
    std::array<Matrix, RANK> products;
    for (int term_index = 0; term_index < RANK; ++term_index) {{
        const Term& term = GERMS[phenotype][term_index];
        Matrix left = linear_form(term.u, ab, h);
        Matrix right = linear_form(term.v, bb, h);
        const std::uint64_t child = splitmix64(node_id ^ (0x9e3779b97f4a7c15ULL * (term_index + 1)));
        products[term_index] = polymorphic_multiply(left, right, h, counter, mode, universe, child);
    }}
    std::array<Matrix, 4> out = {{zero_matrix(h), zero_matrix(h), zero_matrix(h), zero_matrix(h)}};
    for (int term_index = 0; term_index < RANK; ++term_index) {{
        const Term& term = GERMS[phenotype][term_index];
        for (int block = 0; block < 4; ++block) {{
            if (term.w[block] == 0) continue;
            for (std::size_t i = 0; i < out[block].size(); ++i)
                out[block][i] += static_cast<Scalar>(term.w[block]) * products[term_index][i];
        }}
    }}
    return join(out, h);
}}

static Matrix schoolbook(const Matrix& a, const Matrix& b, int n) {{
    Matrix result = zero_matrix(n);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < n; ++j)
                at(result, n, i, j) += at(a, n, i, k) * at(b, n, k, j);
    return result;
}}

static std::uint64_t integer_power(std::uint64_t base, int exponent) {{
    std::uint64_t result = 1;
    for (int i = 0; i < exponent; ++i) result *= base;
    return result;
}}

struct Coverage {{ std::uint64_t covered = 0; std::uint64_t total = 0; }};

static Coverage coverage(const std::vector<int>& selected, const std::vector<int>& universe, int order) {{
    Coverage result;
    const int n = static_cast<int>(universe.size());
    if (order == 1) {{
        for (int i = 0; i < n; ++i) {{
            ++result.total;
            if (any_survives(selected, std::vector<int>{{universe[i]}})) ++result.covered;
        }}
    }} else if (order == 2) {{
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {{
            ++result.total;
            if (any_survives(selected, std::vector<int>{{universe[i], universe[j]}})) ++result.covered;
        }}
    }} else if (order == 3) {{
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
            for (int k = j + 1; k < n; ++k) {{
                ++result.total;
                if (any_survives(selected, std::vector<int>{{universe[i], universe[j], universe[k]}})) ++result.covered;
            }}
    }} else throw std::invalid_argument("unsupported order");
    return result;
}}

int main() {{
    for (int phenotype = 0; phenotype < POOL; ++phenotype)
        if (!verify_tensor_identity(phenotype)) {{
            std::cerr << "invalid phenotype " << phenotype << "\\n";
            return 1;
        }}

    // Exact coefficient-level deduplication.
    std::set<std::array<int, RANK * 12>> distinct;
    for (int p = 0; p < POOL; ++p) {{
        std::array<int, RANK * 12> flat{{}};
        int cursor = 0;
        for (const Term& term : GERMS[p]) {{
            for (int x : term.u) flat[cursor++] = x;
            for (int x : term.v) flat[cursor++] = x;
            for (int x : term.w) flat[cursor++] = x;
        }}
        distinct.insert(flat);
    }}
    if (static_cast<int>(distinct.size()) != POOL) return 2;

    const auto universe = resource_universe();
    if (universe.size() != 44) {{
        std::cerr << "unexpected resource universe size " << universe.size() << "\\n";
        return 3;
    }}
    const std::vector<int> all = [] {{ std::vector<int> x; for (int i = 0; i < POOL; ++i) x.push_back(i); return x; }}();
    const std::vector<int> chosen(SELECTED.begin(), SELECTED.end());
    const Coverage all1 = coverage(all, universe, 1);
    const Coverage all2 = coverage(all, universe, 2);
    const Coverage all3 = coverage(all, universe, 3);
    const Coverage chosen1 = coverage(chosen, universe, 1);
    const Coverage chosen2 = coverage(chosen, universe, 2);
    const Coverage chosen3 = coverage(chosen, universe, 3);
    if (all1.covered != 44 || all1.total != 44 || all2.covered != 871 || all2.total != 946 ||
        all3.covered != 9695 || all3.total != 13244) return 4;
    if (chosen1.covered != 44 || chosen1.total != 44 || chosen2.covered != 622 || chosen2.total != 946 ||
        chosen3.covered != 4874 || chosen3.total != 13244) return 5;

    // Prove the 3-germ cover is cardinality-minimal: no singleton or pair covers all 44 single faults.
    for (int i = 0; i < POOL; ++i)
        if (coverage(std::vector<int>{{i}}, universe, 1).covered == 44) return 6;
    for (int i = 0; i < POOL; ++i) for (int j = i + 1; j < POOL; ++j)
        if (coverage(std::vector<int>{{i, j}}, universe, 1).covered == 44) return 7;

    std::mt19937_64 rng(20260830ULL);
    std::uniform_int_distribution<int> distribution(-2, 2);
    for (const auto& [n, mode] : std::vector<std::pair<int, Mode>>{{{{64, Mode::Mixed}}, {{32, Mode::FaultAware}}}}) {{
        Matrix a = zero_matrix(n), b = zero_matrix(n);
        for (Scalar& x : a) x = distribution(rng);
        for (Scalar& x : b) x = distribution(rng);
        Counter counter;
        const Matrix actual = polymorphic_multiply(a, b, n, counter, mode, universe);
        const Matrix expected = schoolbook(a, b, n);
        int depth = 0; for (int x = n; x > 1; x /= 2) ++depth;
        if (actual != expected || counter.scalar_multiplications != integer_power(7, depth)) return 8;
        std::cout << (mode == Mode::Mixed ? "mixed" : "fault-aware")
                  << " n=" << n << " exact=true scalar_multiplications="
                  << counter.scalar_multiplications << "\\n";
    }}

    std::cout << "phenotypes=" << POOL << " resources=" << universe.size()
              << " minimum_single_fault_cover=" << COVER << "\\n";
    std::cout << "CXX_POOL_INDEPENDENT_VERIFICATION=PASS\\n";
    return 0;
}}
'''
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(source, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    emit(args.certificate, args.output)


if __name__ == "__main__":
    main()
