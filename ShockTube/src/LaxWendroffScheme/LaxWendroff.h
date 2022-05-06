#ifndef LW_H_
#define LW_H_

#include <string>
#include <fstream>
#include "../Riemann.h"

using std::ofstream;
using std::string;

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
        *dMinus, *mMinus, *EMinus;

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

    void Rhs(double F[3], double di, double mi, double Ei)
    {
        double ui = mi / di,
               pi = (Ei - 0.5 * mi * ui) * (gamma - 1.0);

        F[0] = mi;
        F[1] = mi * ui + pi;
        F[2] = ui * (Ei + pi);
    }

    // two-steps Lax-Wendroff scheme
    void UpdateLW(double tx)
    {
        double FPlus[3], F[3], FMinus[3];

        // TODO: 可优化
        for (int i = 1; i < ncell - 1; i++)
        {
            Rhs(FPlus, d[i + 1], m[i + 1], E[i + 1]);
            Rhs(F, d[i], m[i], E[i]);
            Rhs(FMinus, d[i - 1], m[i - 1], E[i - 1]);

            dPlus[i] = 0.5 * (d[i] + d[i + 1]) - 0.5 * tx * (FPlus[0] - F[0]);
            mPlus[i] = 0.5 * (m[i] + m[i + 1]) - 0.5 * tx * (FPlus[1] - F[1]);
            EPlus[i] = 0.5 * (E[i] + E[i + 1]) - 0.5 * tx * (FPlus[2] - F[2]);
            dMinus[i] = 0.5 * (d[i] + d[i - 1]) - 0.5 * tx * (F[0] - FMinus[0]);
            mMinus[i] = 0.5 * (m[i] + m[i - 1]) - 0.5 * tx * (F[1] - FMinus[1]);
            EMinus[i] = 0.5 * (E[i] + E[i - 1]) - 0.5 * tx * (F[2] - FMinus[2]);
        }
    }

    /* time advancing */
    void AdvancingLW(double dt)
    {
        double tx = dt / dx;
        UpdateLW(tx);

        double FPlus[3], FMinus[3];
        for (int i = 1; i < ncell - 1; i++)
        {
            Rhs(FPlus, dPlus[i], mPlus[i], EPlus[i]);
            Rhs(FMinus, dMinus[i], mMinus[i], EMinus[i]);

            d[i] -= tx * (FPlus[0] - FMinus[0]);
            m[i] -= tx * (FPlus[1] - FMinus[1]);
            E[i] -= tx * (FPlus[2] - FMinus[2]);
        }
    }
};

#endif