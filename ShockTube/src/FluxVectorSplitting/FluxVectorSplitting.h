#ifndef FVS_H_
#define FVS_H_

#include "../Riemann.h"

/*
    Piemann solver using flux vector splitting methods (FVS)
    lzmo, ustc, 2022-05-04
*/

// 1D uniform grid for Euler equation
struct Grid
{
    /* data */
    int ncell;
    double xmin, xmax, dx;

    double *x;
    double *d, *m, *E; // conservative variables
    double *dPlus, *mPlus, *EPlus,
        *dMinus, *mMinus, *EMinus; // splitted flux

    Grid(int ncell, double xmin, double xmax)
    {
        this->ncell = ncell;
        this->xmin = xmin;
        this->xmax = xmax;
        dx = (xmax - xmin) / ncell;

        x = new double[ncell];
        d = new double[ncell];
        m = new double[ncell];
        E = new double[ncell];
        dPlus = new double[ncell];
        mPlus = new double[ncell];
        EPlus = new double[ncell];
        dMinus = new double[ncell];
        mMinus = new double[ncell];
        EMinus = new double[ncell];

        for (int i = 0; i < ncell; i++)
        {
            x[i] = dx * (i + 0.5) + xmin;
            d[i] = 0.0;
            m[i] = 0.0;
            E[i] = 0.0;
            dPlus[i] = 0.0;
            mPlus[i] = 0.0;
            EPlus[i] = 0.0;
            dMinus[i] = 0.0;
            mMinus[i] = 0.0;
            EMinus[i] = 0.0;
        }
    }

    ~Grid()
    {
        delete[] x;
        delete[] d;
        delete[] m;
        delete[] E;
        delete[] dPlus;
        delete[] mPlus;
        delete[] EPlus;
        delete[] dMinus;
        delete[] mMinus;
        delete[] EMinus;
    }
};

struct FluxVectorSplitting
{
    /* data */
    Riemann r;

    Grid grid;
};

#endif