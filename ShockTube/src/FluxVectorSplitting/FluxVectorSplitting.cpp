#include "FluxVectorSplitting.h"

bool FluxVectorSplitting::Init(int ncell, double xmin, double xmax, double dL, double uL, double pL, double dR, double uR, double pR)
{
    grid.InitCoordinate(ncell, xmin, xmax);
    grid.InitValue(dL, uL, pL, dR, uR, pR);

    return r.Init(dL, uL, pL, dR, uR, pR);
}

void FluxVectorSplitting::Write(string fileName)
{
    grid.Write(fileName);
}

void FluxVectorSplitting::SplittingSW()
{
    grid.SplittingSW();
}

void FluxVectorSplitting::SplittingVL()
{
    grid.SplittingVL();
}

void FluxVectorSplitting::SplittingLS()
{
    grid.SplittingLS();
}

void FluxVectorSplitting::AdvancingSW(double dt, double tmax)
{
    int step = tmax / dt;

    for (int i = 0; i < step; i++)
    {
        grid.SplittingSW();
        grid.Advancing(dt);
    }
}

void FluxVectorSplitting::AdvancingVL(double dt, double tmax)
{
    int step = tmax / dt;

    for (int i = 0; i < step; i++)
    {
        grid.SplittingVL();
        grid.Advancing(dt);
    }
}

void FluxVectorSplitting::AdvancingLS(double dt, double tmax)
{
    int step = tmax / dt;

    for (int i = 0; i < step; i++)
    {
        grid.SplittingLS();
        grid.Advancing(dt);
    }
}