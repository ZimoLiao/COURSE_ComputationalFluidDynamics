#ifndef POISSON_H_
#define POISSON_H_

#include <iostream>
#include <time.h>
#include "Data.h"

using std::cout;
using std::endl;

// constant value in each band (and P=1)
void TDMA(int number, double W, double *P, double E, double *b, double *x, double &error)
{
    error = 0.0;

    // elimination
    P[0] = 1.0;
    for (int i = 1; i < number; i++)
    {
        P[i] = 1.0 - W * E / P[i - 1];
        b[i] = b[i] - W * b[i - 1] / P[i - 1];
    }

    // back substitution
    x[number - 1] = b[number - 1] / P[number - 1];
    for (int i = number - 2; i >= 0; i--)
    {
        double xn = (b[i] - E * x[i + 1]) / P[i];
        error += pow(xn - x[i], 2.0);
        x[i] = xn;
    }
}

void PoissonJ(Grid &grid, Rhs &rhs, int stepmax, double error)
{
    int num = grid.num, size = grid.size;

    // local data
    double *data = new double[size], *datan = new double[size];
    for (int i = 0; i < size; i++)
    {
        data[i] = grid.data[i];
        datan[i] = grid.data[i];
    }

    // iteration
    clock_t start = clock(), end;
    int step = 0;
    double delta = 1.0;

    while (delta > error && step < stepmax)
    {
        step++;
        delta = 0.0;

        for (int j = 1; j < num; j++)
        {
            for (int i = 1; i < num; i++)
            {
                int indp = Index(i, j, num);
                int inds = Index(i, j - 1, num);
                int indw = Index(i - 1, j, num);
                int indn = Index(i, j + 1, num);
                int inde = Index(i + 1, j, num);

                double phi = (-rhs.data[indp] + data[inds] + data[indw] + data[indn] + data[inde]) / 4.0;

                delta += pow(phi - data[indp], 2.0);
                datan[indp] = phi;
            }
        }

        for (int i = 0; i < size; i++)
        {
            data[i] = datan[i];
        }
    }

    end = clock();
    cout << "time: " << double(end - start) / CLOCKS_PER_SEC << '\n'
         << "delta = " << delta << '\t' << "step = " << step << '\n';

    // copy data
    for (int i = 0; i < size; i++)
    {
        grid.data[i] = data[i];
    }

    // deallocation
    delete[] data;
    delete[] datan;
}

void PoissonGS(Grid &grid, Rhs &rhs, int stepmax, double error)
{
    int num = grid.num, size = grid.size;

    // local data
    double *data = new double[size];
    for (int i = 0; i < size; i++)
    {
        data[i] = grid.data[i];
    }

    // iteration
    clock_t start = clock(), end;
    int step = 0;
    double delta = 1.0;

    while (delta > error && step < stepmax)
    {
        step++;
        delta = 0.0;

        for (int j = 1; j < num; j++)
        {
            for (int i = 1; i < num; i++)
            {
                int indp = Index(i, j, num);
                int inds = Index(i, j - 1, num);
                int indw = Index(i - 1, j, num);
                int indn = Index(i, j + 1, num);
                int inde = Index(i + 1, j, num);

                double phi = (-rhs.data[indp] + data[inds] + data[indw] + data[indn] + data[inde]) / 4.0;

                delta += pow(phi - data[indp], 2.0);
                data[indp] = phi;
            }
        }
    }

    end = clock();
    cout << "time: " << double(end - start) / CLOCKS_PER_SEC << '\n'
         << "delta = " << delta << '\t' << "step = " << step << '\n';

    // copy data
    for (int i = 0; i < size; i++)
    {
        grid.data[i] = data[i];
    }

    // deallocation
    delete[] data;
}

void PoissonSOR(Grid &grid, Rhs &rhs, int stepmax, double error, double omega)
{
    int num = grid.num, size = grid.size;
    double omegac = 1.0 - omega;

    // local data
    double *data = new double[size];
    for (int i = 0; i < size; i++)
    {
        data[i] = grid.data[i];
    }

    // iteration
    clock_t start = clock(), end;
    int step = 0;
    double delta = 1.0;

    while (delta > error && step < stepmax)
    {
        step++;
        delta = 0.0;

        for (int j = 1; j < num; j++)
        {
            for (int i = 1; i < num; i++)
            {
                int indp = Index(i, j, num);
                int inds = Index(i, j - 1, num);
                int indw = Index(i - 1, j, num);
                int indn = Index(i, j + 1, num);
                int inde = Index(i + 1, j, num);

                double phi = (-rhs.data[indp] + data[inds] + data[indw] + data[indn] + data[inde]) / 4.0;

                delta += pow(phi - data[indp], 2.0);

                // update
                data[indp] *= omegac;
                data[indp] += omega * phi;
            }
        }
    }

    end = clock();
    cout << "time: " << double(end - start) / CLOCKS_PER_SEC << '\n'
         << "delta = " << delta << '\t' << "step = " << step << '\n';

    // copy data
    for (int i = 0; i < size; i++)
    {
        grid.data[i] = data[i];
    }

    // deallocation
    delete[] data;
}

// Line G-S method (in Y-direction, i.e. colume)
void PoissonLGSY(Grid &grid, Rhs &rhs, int stepmax, double error)
{
    int num = grid.num, size = grid.size;

    // local data
    double *data = new double[size];
    for (int i = 0; i < size; i++)
    {
        data[i] = grid.data[i];
    }

    // iteration
    clock_t start = clock(), end;
    int step = 0;
    double delta = 1.0;

    int number = num - 1; // number of points in each direction
    double W = -0.25, E = -0.25;
    double *P = new double[number], *b = new double[number];
    double delta_tdma;

    while (delta > error && step < stepmax)
    {
        step++;
        delta = 0.0;

        for (int j = 1; j < num; j++)
        {
            int ind = Index(1, j, num);

            for (int i = 0; i < number; i++)
            {
                b[i] = 0.25 * (-rhs.data[ind + i] + data[Index(i + 1, j - 1, num)] + data[Index(i + 1, j + 1, num)]);
            }

            TDMA(number, W, P, E, b, &data[ind], delta_tdma);

            delta += delta_tdma;
        }
    }

    end = clock();
    cout << "time: " << double(end - start) / CLOCKS_PER_SEC << '\n'
         << "delta = " << delta << '\t' << "step = " << step << '\n';

    // copy data
    for (int i = 0; i < size; i++)
    {
        grid.data[i] = data[i];
    }

    // deallocation
    delete[] P;
    delete[] b;
    delete[] data;
}

void PoissonMG_GS(Grid &grid, Rhs &rhs, int n1, int n2, double error)
{
    int num = grid.num, numc = num / 2;
    cout << num << endl; // TODO

    if (num % 2 || num < 50)
    {
        PoissonGS(grid, rhs, 10000, error);
    }
    else
    {
        PoissonGS(grid, rhs, n1, error);

        rhs.Residual(grid);
        Rhs rhsf(rhs), rhsc(numc);
        rhsc.MGrestriction(rhs);

        Grid gridc(numc);

        PoissonMG_GS(gridc, rhsc, n1, n2, error / 10);

        cout << num << endl; // TODO
        Grid grid_interp(num);
        grid_interp.MGinterpolation(gridc);

        PoissonGS(grid_interp, rhsf, n2, error);

        grid += grid_interp;
    }
}

#endif