#include "Grid.h"

constexpr int NCELL = 128;

int main()
{
    Info info;
    info.error = 1e-6;

    int level = 4;
    double omega = 1.5; // TODO: omega的选择！
    // omega = 2 / (1 + sin(3.1415926535897932 / NCELL));

    Grid phi(NCELL), rhs, res(NCELL);

    rhs.InitPoissonRhs(NCELL);

    /* Jacobi method */
    {
        phi.Flush();
        res.Flush();
        info.Reset();

        phi.SolvePoissonJacobi(rhs, info);
        phi.Write("Jacobi-sol.plt", "phi");

        res.InitPoissonRes(phi, rhs);
        res.Write("Jacobi-res.plt", "phi_res");

        info.Print("Jacobi");
        info.Write("Jacobi-info.plt");
    }

    /* G-S method */
    {
        phi.Flush();
        res.Flush();
        info.Reset();

        phi.SolvePoissonGS(rhs, info);
        phi.Write("GS-sol.plt", "phi");

        res.InitPoissonRes(phi, rhs);
        res.Write("GS-res.plt", "phi_res");

        info.Print("GS");
        info.Write("GS-info.plt");
    }

    /* SOR method */
    {
        phi.Flush();
        res.Flush();
        info.Reset();

        phi.SolvePoissonSOR(rhs, omega, info);
        // phi.SolvePoissonSOR(rhs, info);
        phi.Write("SOR-sol.plt", "phi");

        res.InitPoissonRes(phi, rhs);
        res.Write("SOR-res.plt", "phi_res");

        info.Print("SOR");
        info.Write("SOR-info.plt");
    }

    /* Line G-S method */
    {
        phi.Flush();
        res.Flush();
        info.Reset();

        phi.SolvePoissonLGSY(rhs, info);
        phi.Write("LGSY-sol.plt", "phi");

        res.InitPoissonRes(phi, rhs);
        res.Write("LGSY-res.plt", "phi_res");

        info.Print("LGSY");
        info.Write("LGSY-info.plt");
    }

    /* Multi-grid G-S method */
    {
        phi.Flush();
        res.Flush();
        info.Reset();

        phi.SolvePoissonMGGS(rhs, level, 2, 2, info);
        phi.Write("MG-GS-sol.plt", "phi");

        res.InitPoissonRes(phi, rhs);
        res.Write("MG-GS-res.plt", "phi_res");

        info.Print("MG-GS");
        info.Write("MG-GS-info.plt");
    }

    /* Multi-grid SOR method */
    {
        phi.Flush();
        res.Flush();
        info.Reset();

        phi.SolvePoissonMGSOR(rhs, omega, level, 2, 2, info);
        // phi.SolvePoissonMGSOR(rhs, level, 2, 2, info);
        phi.Write("MG-SOR-sol.plt", "phi");

        res.InitPoissonRes(phi, rhs);
        res.Write("MG-SOR-res.plt", "phi_res");

        info.Print("MG-SOR");
        info.Write("MG-SOR-info.plt");
    }
}