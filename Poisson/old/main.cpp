#include "Poisson.h"

constexpr double pi = 3.14159265358979;

int main()
{
    int num = 100;
    double omega = 1.0;

    Grid jacobi(num), gs(num), sor(num), lgs(num);
    Rhs rhs(num);

    // PoissonJ(jacobi, rhs, 10000, 1e-6);
    // jacobi.WriteAscii("Jacobi.plt");

    clock_t start = clock(), end;
    PoissonGS(gs, rhs, 10000, 1e-6);
    end = clock();
    double gstime = double(end - start) / CLOCKS_PER_SEC;
    gs.WriteAscii("GS.plt");

    // omega = 2 / (1 + sin(pi / num));
    // PoissonSOR(sor, rhs, 10000, 1e-6, omega);
    // sor.WriteAscii("SOR.plt");

    // PoissonLGSY(lgs, rhs, 10000, 1e-6);
    // lgs.WriteAscii("LGSY.plt");

    /* Multi-grid test */
    cout << "\nMulti-grid test\n";

    // TODO: test 1
    // Grid coarse(num / 2), fine(num);

    // coarse.MGrestriction(sor);
    // coarse.WriteAscii("SOR_coarse.plt");

    // fine.MGinterpolation(coarse);
    // fine.WriteAscii("SOR_fine.plt");

    // TODO: test 2
    // Grid fine(num), coarse(num / 2);
    // PoissonGS(fine, rhs, 500, 1e-4);
    // fine.WriteAscii("MG-GS-approximation.plt");

    // rhs.Residual(fine);
    // Rhs rhsf(rhs);
    // Rhs rhs_coarse(num / 2);
    // rhs_coarse.MGrestriction(rhs);

    // PoissonGS(coarse, rhs_coarse, 500, 1e-8);

    // Grid rec(num);
    // rec.MGinterpolation(coarse);

    // PoissonGS(rec, rhsf, 500, 1e-6);
    // rec.WriteAscii("correction.plt");

    // fine += rec;
    // fine.WriteAscii("MG-GS-solution.plt");

    Grid mggs(num);

    start = clock();
    cout << "cycle\n\n";
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    cout << "cycle\n\n";
    rhs.Init(num);
    PoissonMG_GS(mggs, rhs, 10, 10, 1e-6);
    end = clock();
    double mggstime = double(end - start) / CLOCKS_PER_SEC;

    mggs.WriteAscii("MG-GS.plt");

    cout << "\nGS time: " << gstime << "\nMG-GS time: " << mggstime << '\n';
}