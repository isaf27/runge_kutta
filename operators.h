#pragma once

#include <array>
#include <ostream>

template <size_t N, typename Real>
std::ostream& operator<<(std::ostream& out, const std::array<Real, N>& vect) {
    if (N == 1) {
        return out << vect[0];
    }
    out << "(";
    for (size_t i = 0; i < N; ++i) {
        if (i > 0) {
            out << ", ";
        }
        out << vect[i];
    }
    out << ")";
    return out;
}

template <size_t N, typename Real>
std::array<Real, N> operator+(const std::array<Real, N>& first, const std::array<Real, N>& second) {
    std::array<Real, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = first[i] + second[i];
    }
    return result;
}

template <size_t N, typename Real, typename Coef>
std::array<Real, N> operator*(const std::array<Real, N>& vect, const Coef& coef) {
    std::array<Real, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = vect[i] * coef;
    }
    return result;
}