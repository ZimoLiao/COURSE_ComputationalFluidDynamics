#include <iostream>
#include "Numerical.h"
#include "Eno.h"
#include "Weno.h"

using namespace std;

int main()
{
    // ENO approximation test
    /*
    constexpr int ncell = 10;
    Grid x(ncell, 'm'), v_mean(ncell, 'c'), p(ncell, 'm');
    x.InitUniformCoordinate(-1, 1);
    v_mean(1) = 1;
    v_mean(2) = 2;
    v_mean(3) = 3;
    v_mean(4) = 4;
    v_mean(5) = 5;
    v_mean(6) = -1;
    v_mean(7) = -1;
    v_mean(8) = -1;
    v_mean(9) = -1;
    v_mean(10) = -1;

    Eno(x, v_mean, p, 4);

    cout << "u(x) - cell mean value:\n";
    v_mean.Print();
    cout << "v(x+-) -cell boundary approximation (ENO):\n";
    p.Print();
    */

    // Lagrangian interpolation test
    /*
    double x[3] = {1, 2, 3}, xp = 1.5, c[3];
    InterpL(x, 3, xp, c);
    */

    // WENO approximation
    constexpr int ncell = 10;
    Grid x(ncell, 'm'), v_mean(ncell, 'c'), p(ncell, 'f');
    x.InitUniformCoordinate(-1, 1);
    v_mean(1) = 1;
    v_mean(2) = 4;
    v_mean(3) = 9;
    v_mean(4) = 16;
    v_mean(5) = 25;
    v_mean(6) = 36;
    v_mean(7) = 0;
    v_mean(8) = 0;
    v_mean(9) = 0;
    v_mean(10) = 0;
    Weno_JS(x, v_mean, p, 2);
    v_mean.Print();
    p.Print();
}