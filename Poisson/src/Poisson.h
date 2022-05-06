#ifndef POISSON_H_
#define POISSON_H_

#include <iostream>
#include <time.h>
#include "Grid.h"

using std::cout;
using std::endl;

constexpr int STEPMAX = 100000;

/* information of iterative methods */
struct Info
{
    /* data */
    double time = 0.0, error = 1e-6;
    int step = 0;

    double iterErrorN2[STEPMAX], iterErrorNinf[STEPMAX];

    /* functions */
    void Start()
    {
        start = clock();
    }

    void End()
    {
        end = clock();
        time = double(end - start) / CLOCKS_PER_SEC;
    }

private:
    clock_t start, end;
};

/* iterative solvers */
void PoissonGS(Grid &phi, Grid &rhs, double error, Info &info)
{
    int nCell = phi.nCell, nSize = phi.nSize;

    // local data
    double *data = new double[nSize];
    for (int ind = 0; ind < nSize; ind++)
    {
        data[ind] = phi.data[ind];
    }

    // iteration
    info.Start();

    double error = 1.0;
    while (info.step < STEPMAX && error > info.error)
    {
        info.step++;
        error = 0.0;

        for (int j = 1; j < nCell; j++)
        {
            for (int i = 1; i < nCell; i++)
            {
            }
        }
    }

    info.End();

    delete[] phi;
}

#endif