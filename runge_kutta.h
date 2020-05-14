#pragma once

#include <functional>
#include <vector>
#include <algorithm>

template <typename Real, typename Vect>
std::vector<Vect> runge_kutta_1(
    Real x0,
    Vect y0,
    std::function<Vect(Real, Vect)> f,
    Real step,
    size_t l,
    size_t r
) {
    if (l > r) {
        return std::vector<Vect>();
    }
    std::vector<Vect> result(r - l + 1);
    Real x = x0;
    Vect y = y0;
    if (l == 0) {
        result[0] = y;
    }
    for (size_t i = 0; i < r; ++i) {
        Vect k1 = f(x, y);
        x = x + step;
        y = y + k1 * step;
        if (l <= i + 1) {
            result[i + 1 - l] = y;
        }
    }
    return result;
}

template <typename Real, typename Vect>
std::vector<Vect> runge_kutta_2(
    Real x0,
    Vect y0,
    std::function<Vect(Real, Vect)> f,
    Real step,
    size_t l,
    size_t r
) {
    if (l > r) {
        return std::vector<Vect>();
    }
    std::vector<Vect> result(r - l + 1);
    Real x = x0;
    Vect y = y0;
    if (l == 0) {
        result[0] = y;
    }
    for (size_t i = 0; i < r; ++i) {
        Vect k1 = f(x, y);
        Vect k2 = f(x + step, y + k1 * step);
        x = x + step;
        y = y + (k1 + k2) * (step / (Real)2);
        if (l <= i + 1) {
            result[i + 1 - l] = y;
        }
    }
    return result;
}

template <typename Real, typename Vect>
std::vector<Vect> runge_kutta_3(
    Real x0,
    Vect y0,
    std::function<Vect(Real, Vect)> f,
    Real step,
    size_t l,
    size_t r
) {
    if (l > r) {
        return std::vector<Vect>();
    }
    std::vector<Vect> result(r - l + 1);
    Real x = x0;
    Vect y = y0;
    if (l == 0) {
        result[0] = y;
    }
    for (size_t i = 0; i < r; ++i) {
        Vect k1 = f(x, y);
        Vect k2 = f(x + (step / (Real)2), y + k1 * (step / (Real)2));
        Vect k3 = f(x + step, y + k1 * (-step) + k2 * ((Real)2 * step));
        x = x + step;
        y = y + (k1 + k2 * (Real)4 + k3) * (step / (Real)6);
        if (l <= i + 1) {
            result[i + 1 - l] = y;
        }
    }
    return result;
}

template <typename Real, typename Vect>
std::vector<Vect> runge_kutta_4(
    Real x0,
    Vect y0,
    std::function<Vect(Real, Vect)> f,
    Real step,
    size_t l,
    size_t r
) {
    if (l > r) {
        return std::vector<Vect>();
    }
    std::vector<Vect> result(r - l + 1);
    Real x = x0;
    Vect y = y0;
    if (l == 0) {
        result[0] = y;
    }
    for (size_t i = 0; i < r; ++i) {
        Vect k1 = f(x, y);
        Vect k2 = f(x + (step / (Real)2), y + k1 * (step / (Real)2));
        Vect k3 = f(x + (step / (Real)2), y + k2 * (step / (Real)2));
        Vect k4 = f(x + step, y + k3 * step);
        x = x + step;
        y = y + (k1 + k2 * (Real)2 + k3 * (Real)2 + k4) * (step / (Real)6);
        if (l <= i + 1) {
            result[i + 1 - l] = y;
        }
    }
    return result;
}