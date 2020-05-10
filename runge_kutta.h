#pragma once

#include <functional>
#include <vector>
#include <algorithm>

template <typename Real, typename Vect>
Vect runge_kutta(
    Real x0,
    Vect y0,
    std::function<Vect(Real, Vect)> f,
    Real step,
    size_t iterations
) {
    Real x = x0;
    Vect y = y0;
    for (size_t i = 0; i < iterations; ++i) {
        Vect k1 = f(x, y);
        Vect k2 = f(x + (step / (Real)2), y + k1 * (step / (Real)2));
        Vect k3 = f(x + (step / (Real)2), y + k2 * (step / (Real)2));
        Vect k4 = f(x + step, y + k3 * step);
        x = x + step;
        y = y + (k1 + k2 * (Real)2 + k3 * (Real)2 + k4) * (step / (Real)6);
    }
    return y;
}

template <typename Real, typename Vect>
std::vector<Vect> runge_kutta_range(
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

template <typename Real, typename Vect>
Vect runge_kutta(Real x0, Vect y0, Real x_need, std::function<Vect(Real, Vect)> f, Real step) {
    if ((step > 0) ^ (x_need > x0)) {
        step = -step;
    }
    size_t iterations = (x_need - x0) / step;
    Real x_l = x0 + step * iterations;
    Real x_r = x0 + step * (iterations + 1);
    std::vector<Vect> val = runge_kutta_range(x0, y0, f, step, iterations, iterations + 1);
    return val[0] * ((x_r - x_need) / (x_r - x_l)) + val[1] * ((x_need - x_l) / (x_r - x_l));
}