#include "operators.h"
#include "runge_kutta.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using ld = long double;
using Vect = std::array<ld, 4>;

constexpr ld c3 = (ld)1 / (ld)5;
constexpr ld c2 = (ld)7401 / (ld)100;
constexpr ld c1 = 5;
constexpr ld c0 = (ld)4901 / (ld)4;

Vect f(ld x, Vect y) {
    Vect result;
    result[0] = y[1];
    result[1] = y[2];
    result[2] = y[3];
    result[3] = -c3 * y[3] - c2 * y[2] - c1 * y[1] - c0 * y[0];
    return result;
}

ld value(ld x) {
    return 5 * exp(-x / (ld)10) * cos(7 * x) + 3 * sin(5 * x);
}

constexpr ld X0 = 0;
constexpr Vect Y0 = Vect({5, (ld)29 / (ld)2, -(ld)4899 / (ld)20, -(ld)60301 / (ld)200});
constexpr ld STEP = 0.05;
constexpr size_t ITER = 400;

int main(int argc, char *argv[]) {
    std::ofstream fout("report2.csv");
    fout << std::fixed << std::setprecision(10);
    std::vector<Vect> res1 = runge_kutta<ld, Vect>(X0, Y0, f, STEP, 0, ITER, Method::One);
    std::vector<Vect> res2 = runge_kutta<ld, Vect>(X0, Y0, f, STEP, 0, ITER, Method::Two);
    std::vector<Vect> res3 = runge_kutta<ld, Vect>(X0, Y0, f, STEP, 0, ITER, Method::Three);
    std::vector<Vect> res4 = runge_kutta<ld, Vect>(X0, Y0, f, STEP, 0, ITER, Method::Four);
    fout << "x,y,y (p = 1),error (p = 1),y (p = 2),error (p = 2),y (p = 3),error (p = 3),y (p = 4),error (p = 4)\n";
    for (size_t i = 0; i <= ITER; i++) {
        ld real = value(X0 + i * STEP);
        fout << X0 + i * STEP << ",";
        fout << real << ",";
        fout << res1[i][0] << ",";
        fout << res1[i][0] - real << ",";
        fout << res2[i][0] << ",";
        fout << res2[i][0] - real << ",";
        fout << res3[i][0] << ",";
        fout << res3[i][0] - real << ",";
        fout << res4[i][0] << ",";
        fout << res4[i][0] - real << "\n";
    }
    fout.close();
    return 0;
}