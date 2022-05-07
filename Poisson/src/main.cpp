#include "Grid.h"

constexpr int NCELL = 128;

int main()
{
    Info info;
    info.error = 1e-6;

    Grid phi(NCELL), rhs, res(NCELL);

    rhs.InitPoissonRhs(NCELL);

    /* G-S method */
    phi.Flush();
    res.Flush();
    info.Reset();

    phi.SolvePoissonGS(rhs, info);
    phi.Write("GS-sol.plt", "phi");

    res.InitPoissonRes(phi, rhs);
    res.Write("GS-res.plt", "phi_res");

    info.Print("GS");
    info.Write("GS-info.plt");

    /* SOR method */
    phi.Flush();
    res.Flush();
    info.Reset();

    double omega = 1.5;
    phi.SolvePoissonSOR(rhs, omega, info);
    phi.Write("SOR-sol.plt", "phi");

    res.InitPoissonRes(phi, rhs);
    res.Write("SOR-res.plt", "phi_res");

    info.Print("SOR");
    info.Write("SOR-info.plt");

    /* Multi-grid G-S method */
    phi.Flush();
    res.Flush();
    info.Reset();

    phi.SolvePoissonMGGS(rhs, 4, 2, 2, info);
    phi.Write("MG-GS-sol.plt", "phi");

    res.InitPoissonRes(phi, rhs);
    res.Write("MG-GS-res.plt", "phi_res");

    info.Print("MG-GS");
    info.Write("MG-GS-info.plt");

    /* Multi-grid SOR method */
    phi.Flush();
    res.Flush();
    info.Reset();

    phi.SolvePoissonMGSOR(rhs, omega, 4, 2, 2, info);
    phi.Write("MG-SOR-sol.plt", "phi");

    res.InitPoissonRes(phi, rhs);
    res.Write("MG-SOR-res.plt", "phi_res");

    info.Print("MG-SOR");
    info.Write("MG-SOR-info.plt");
}