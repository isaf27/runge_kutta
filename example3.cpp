#include "runge_kutta.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using ld = long double;

ld f(ld x, ld y) {
    ld result = 0;
    result += 5 * sin(x * x + y) * cos(y * y * y + x * x);
    result += (log(2 * x * x + 3) * (x + y) * (x + y)) / (ld)2;
    result -= (2 * exp(-x / 2) * sqrtl(x * x + y * y * y * y)) / (1 + x * x);
    return result;
}

constexpr ld X0 = 0;
constexpr ld Y0 = 0;
constexpr ld STEP = 0.01;
constexpr size_t ITER = 150;

int main(int argc, char *argv[]) {
    std::ofstream fout("report3.csv");
    fout << std::fixed << std::setprecision(10);
    std::vector<ld> res4 = runge_kutta<ld, ld>(X0, Y0, f, STEP, 0, ITER, Method::Four);
    fout << "x,y\n";
    for (size_t i = 0; i <= ITER; i++) {
        fout << X0 + i * STEP << ",";
        fout << res4[i] << "\n";
    }
    fout.close();
    return 0;
}