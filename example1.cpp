#include "operators.h"
#include "runge_kutta.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using ld = long double;
using Vect = std::array<ld, 2>;

Vect f(ld x, Vect y) {
    Vect result;
    result[0] = y[1];
    result[1] = y[1] - (1 + x) * (1 + tan(y[0])) * y[1] * y[1];
    return result;
}

constexpr ld X0 = 0;
constexpr Vect Y0 = Vect({0, 1});
constexpr ld STEP = 0.2;
constexpr size_t ITER = 25;

int main(int argc, char *argv[]) {
    std::ofstream fout("report1.csv");
    fout << std::fixed << std::setprecision(10);
    std::vector<Vect> res1 = runge_kutta<ld, Vect>(X0, Y0, f, STEP, 0, ITER, Method::One);
    std::vector<Vect> res2 = runge_kutta<ld, Vect>(X0, Y0, f, STEP, 0, ITER, Method::Two);
    std::vector<Vect> res3 = runge_kutta<ld, Vect>(X0, Y0, f, STEP, 0, ITER, Method::Three);
    std::vector<Vect> res4 = runge_kutta<ld, Vect>(X0, Y0, f, STEP, 0, ITER, Method::Four);
    fout << "x,y,y (p = 1),error (p = 1),y (p = 2),error (p = 2),y (p = 3),error (p = 3),y (p = 4),error (p = 4)\n";
    for (size_t i = 0; i <= ITER; ++i) {
        ld real = atanl(X0 + i * STEP);
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