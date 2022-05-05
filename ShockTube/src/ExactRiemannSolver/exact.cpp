#include <fstream>
#include "ExactRiemannSolver.h"

using std::ofstream;

constexpr int n = 100;

int main()
{
    ExactRiemannSolver ers;

    ers.Init(1.0, 0.0, 1.0, 0.125, 0.0, 0.1); // Sod problem
    // ers.Init(1.0, 0.0, 1000.0, 1.0, 0.0, 0.01); // blast wave problem

    double smax = 2.0;
    // double smax = 125.0 / 3.0;

    double s[n], d[n], u[n], p[n];
    double ds = 2 * smax / n;
    for (int i = 0; i < n; i++)
    {
        s[i] = -smax + ds * i;
    }

    ers.Solve(n, s, d, u, p);
    ers.Write("exact.plt", n, s, d, u, p);
}