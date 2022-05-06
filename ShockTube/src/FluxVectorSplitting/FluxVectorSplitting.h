#ifndef FVS_H_
#define FVS_H_

#include <string>
#include <fstream>
#include "../Riemann.h"

using std::ofstream;
using std::string;

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
    double gamma = 1.4;

    double *x;
    double *d, *m, *E; // conservative variables
    double *dPlus, *mPlus, *EPlus,
        *dMinus, *mMinus, *EMinus; // splitted flux

    Grid()
    {
        x = new double[1];
        d = new double[1];
        m = new double[1];
        E = new double[1];
        dPlus = new double[1];
        mPlus = new double[1];
        EPlus = new double[1];
        dMinus = new double[1];
        mMinus = new double[1];
        EMinus = new double[1];
    }

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

    void InitCoordinate(int ncell, double xmin, double xmax)
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

    void InitValue(double dL, double uL, double pL, double dR, double uR, double pR)
    {
        double mL = uL * dL,
               EL = dL * (0.5 * uL * uL + pL / ((gamma - 1.0) * dL)),
               mR = uR * dR,
               ER = dR * (0.5 * uR * uR + pR / ((gamma - 1.0) * dR));

        for (int i = 0; i < ncell; i++)
        {
            if (x[i] < 0.0)
            { // left
                d[i] = dL;
                m[i] = mL;
                E[i] = EL;
            }
            else
            { // right
                d[i] = dR;
                m[i] = mR;
                E[i] = ER;
            }
        }
    }

    void Write(string fileName)
    {
        ofstream fout;
        fout.open(fileName);

        fout << "variables = x d u p e\n"; // including specific internal energy

        for (int i = 0; i < ncell; i++)
        {
            double u = m[i] / d[i],
                   p = (E[i] - 0.5 * m[i] * u) * (gamma - 1.0);
            fout << x[i] << '\t' << d[i] << '\t' << u << '\t' << p << '\t' << p / d[i] / (gamma - 1.0) << '\n';
        }

        fout.close();
    }

    /* splitting schemes */
    // Steger-Warming splitting
    void SplittingSW()
    {
        for (int i = 0; i < ncell; i++)
        {
            double u = m[i] / d[i],
                   p = (E[i] - 0.5 * m[i] * u) * (gamma - 1.0),
                   a = sqrt(gamma * p / d[i]);

            double lambda1 = u - a, lambda2 = u, lambda3 = u + a;

            double lambda1Plus, lambda2Plus, lambda3Plus,
                lambda1Minus, lambda2Minus, lambda3Minus;

            double H = (E[i] + p) / d[i];

            lambda1Plus = 0.5 * (lambda1 + fabs(lambda1));
            lambda2Plus = 0.5 * (lambda2 + fabs(lambda2));
            lambda3Plus = 0.5 * (lambda3 + fabs(lambda3));
            lambda1Minus = 0.5 * (lambda1 - fabs(lambda1));
            lambda2Minus = 0.5 * (lambda2 - fabs(lambda2));
            lambda3Minus = 0.5 * (lambda3 - fabs(lambda3));

            dPlus[i] = lambda1Plus + 2.0 * (gamma - 1.0) * lambda2Plus + lambda3Plus;
            mPlus[i] = lambda1 * lambda1Plus + 2.0 * (gamma - 1.0) * lambda2Plus + lambda3 * lambda3Plus;
            EPlus[i] = (H - u * a) * lambda1Plus + (gamma - 1.0) * u * u * lambda2Plus + (H + u * a) * lambda3Plus;

            dMinus[i] = lambda1Minus + 2.0 * (gamma - 1.0) * lambda2Minus + lambda3Minus;
            mMinus[i] = lambda1 * lambda1Minus + 2.0 * (gamma - 1.0) * lambda2Minus + lambda3 * lambda3Minus;
            EMinus[i] = (H - u * a) * lambda1Minus + (gamma - 1.0) * u * u * lambda2Minus + (H + u * a) * lambda3Minus;

            double coef = d[i] / (2.0 * gamma);
            dPlus[i] *= coef;
            mPlus[i] *= coef;
            EPlus[i] *= coef;
            dMinus[i] *= coef;
            mMinus[i] *= coef;
            EMinus[i] *= coef;
        }
    }

    // van Leer splitting
    void SplittingVL()
    {
        for (int i = 0; i < ncell; i++)
        {
            double u = m[i] / d[i],
                   p = (E[i] - 0.5 * m[i] * u) * (gamma - 1.0),
                   a = sqrt(gamma * p / d[i]),
                   M = u / a;

            double termPlus = 0.5 * (gamma - 1.0) * M + 1.0,
                   termMinus = 0.5 * (gamma - 1.0) * M - 1.0,
                   coefPlus = 0.25 * d[i] * a * pow(1.0 + M, 2.0),
                   coefMinus = -0.25 * d[i] * a * pow(1.0 - M, 2.0);

            dPlus[i] = 1.0;
            mPlus[i] = 2.0 * a / gamma * termPlus;
            EPlus[i] = 2.0 * a * a / (gamma * gamma - 1.0) * termPlus * termPlus;
            dPlus[i] *= coefPlus;
            mPlus[i] *= coefPlus;
            EPlus[i] *= coefPlus;

            dMinus[i] = 1.0;
            mMinus[i] = 2.0 * a / gamma * termMinus;
            EMinus[i] = 2.0 * a * a / (gamma * gamma - 1.0) * termMinus * termMinus;
            dMinus[i] *= coefMinus;
            mMinus[i] *= coefMinus;
            EMinus[i] *= coefMinus;
        }
    }

    // Liou-Steffen scheme
    // also refers to the advection upstream splitting method (AUSM)
    void SplittingLS()
    {
    }

    /* time advancing */
    void Advancing(double dt)
    {
        for (int i = 1; i < ncell - 1; i++)
        {
            d[i] -= dt / dx * ((dPlus[i] + dMinus[i + 1]) - (dPlus[i - 1] + dMinus[i]));
            m[i] -= dt / dx * ((mPlus[i] + mMinus[i + 1]) - (mPlus[i - 1] + mMinus[i]));
            E[i] -= dt / dx * ((EPlus[i] + EMinus[i + 1]) - (EPlus[i - 1] + EMinus[i]));
        }
    }
};

struct FluxVectorSplitting
{
    /* data */
    Riemann r;

    Grid grid;

    /* functions */
    bool Init(int ncell, double xmin, double xmax, double dL, double uL, double pL, double dR, double uR, double pR);
    void Write(string fileName);

    /* splitting schemes */
    // Steger-Warming splitting
    void SplittingSW();

    // van Leer splitting
    void SplittingVL();

    // Liou-Steffen scheme
    // also refers to the advection upstream splitting method (AUSM)
    void SplittingLS();

    /* time advancing */
    void AdvancingSW(double dt, double tmax);
    void AdvancingVL(double dt, double tmax);
    void AdvancingLS(double dt, double tmax);
};

#endif