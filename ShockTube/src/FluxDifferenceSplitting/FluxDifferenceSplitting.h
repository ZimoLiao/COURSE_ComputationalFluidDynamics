#ifndef FDS_H_
#define FDS_H_

#include <string>
#include <fstream>
#include "../Riemann.h"

using std::ofstream;
using std::string;

/*
    Piemann solver using flux difference splitting methods (FDS, or Godunov's method)
    lzmo, ustc, 2022-05-10
*/

// 1D uniform grid for Euler equation
struct Grid
{
    /* data */
    int ncell;
    double xmin, xmax, dx;
    double gamma = 1.4;

    double *x;
    double *d, *m, *E;                 // conservative variables
    double *dFlux, *mFlux, *EFlux;     // flux
    double *dRoe, *uRoe, *hRoe, *aRoe; // Roe average

    Grid()
    {
        x = new double[1];
        d = new double[1];
        m = new double[1];
        E = new double[1];
        dFlux = new double[1];
        mFlux = new double[1];
        EFlux = new double[1];
        dRoe = new double[1];
        uRoe = new double[1];
        hRoe = new double[1];
        aRoe = new double[1];
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
        dFlux = new double[ncell - 1];
        mFlux = new double[ncell - 1];
        EFlux = new double[ncell - 1];
        dRoe = new double[ncell - 1];
        uRoe = new double[ncell - 1];
        hRoe = new double[ncell - 1];
        aRoe = new double[ncell - 1];

        for (int i = 0; i < ncell; i++)
        {
            x[i] = dx * (i + 0.5) + xmin;
            d[i] = 0.0;
            m[i] = 0.0;
            E[i] = 0.0;
            if (i < ncell - 1)
            {
                dFlux[i] = 0.0;
                mFlux[i] = 0.0;
                EFlux[i] = 0.0;
                dRoe[i] = 0.0;
                uRoe[i] = 0.0;
                hRoe[i] = 0.0;
                aRoe[i] = 0.0;
            }
        }
    }

    ~Grid()
    {
        delete[] x;
        delete[] d;
        delete[] m;
        delete[] E;
        delete[] dFlux;
        delete[] mFlux;
        delete[] EFlux;
        delete[] dRoe;
        delete[] uRoe;
        delete[] hRoe;
        delete[] aRoe;
    }

    void InitCoordinate(int ncell, double xmin, double xmax)
    {
        this->ncell = ncell;
        this->xmin = xmin;
        this->xmax = xmax;
        dx = (xmax - xmin) / ncell;

        delete[] x;
        delete[] d;
        delete[] m;
        delete[] E;
        delete[] dFlux;
        delete[] mFlux;
        delete[] EFlux;
        delete[] dRoe;
        delete[] uRoe;
        delete[] hRoe;
        delete[] aRoe;
        x = new double[ncell];
        d = new double[ncell];
        m = new double[ncell];
        E = new double[ncell];
        dFlux = new double[ncell - 1];
        mFlux = new double[ncell - 1];
        EFlux = new double[ncell - 1];
        dRoe = new double[ncell - 1];
        uRoe = new double[ncell - 1];
        hRoe = new double[ncell - 1];
        aRoe = new double[ncell - 1];

        for (int i = 0; i < ncell; i++)
        {
            x[i] = dx * (i + 0.5) + xmin;
            d[i] = 0.0;
            m[i] = 0.0;
            E[i] = 0.0;
            if (i < ncell - 1)
            {
                dFlux[i] = 0.0;
                mFlux[i] = 0.0;
                EFlux[i] = 0.0;
                dRoe[i] = 0.0;
                uRoe[i] = 0.0;
                hRoe[i] = 0.0;
                aRoe[i] = 0.0;
            }
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

    /* flux approximation */
    // Roe-Pike method
    void FluxRP()
    {
        double dLsr, uL, hL, pL, dRsr, uR, hR, pR;
        double dDelta, pDelta, uDelta;
        double lambda1, lambda2, lambda3;
        // int sign1, sign2, sign3;
        double alpha1, alpha2, alpha3;
        double dFluxL, mFluxL, EFluxL;
        double dFluxR, mFluxR, EFluxR;

        for (int i = 0; i < ncell - 1; i++)
        {
            // Roe average value
            dLsr = sqrt(d[i]);
            uL = m[i] / d[i];
            pL = (gamma - 1.0) * (E[i] - 0.5 * m[i] * m[i] / d[i]);
            hL = (E[i] + pL) / d[i];

            dRsr = sqrt(d[i + 1]);
            uR = m[i + 1] / d[i + 1];
            pR = (gamma - 1.0) * (E[i + 1] - 0.5 * m[i + 1] * m[i + 1] / d[i + 1]);
            hR = (E[i + 1] + pR) / d[i + 1];

            dDelta = d[i + 1] - d[i];
            uDelta = uR - uL;
            pDelta = pR - pL;

            dRoe[i] = dLsr * dRsr;
            uRoe[i] = (dLsr * uL + dRsr * uR) / (dLsr + dRsr);
            hRoe[i] = (dLsr * hL + dRsr * hR) / (dLsr + dRsr);
            aRoe[i] = sqrt((gamma - 1.0) * (hRoe[i] - 0.5 * uRoe[i] * uRoe[i]));

            // eigenvalues
            lambda1 = uRoe[i] - aRoe[i];
            lambda2 = uRoe[i];
            lambda3 = uRoe[i] + aRoe[i];

            // sign1 = lambda1 <= 0;
            // // sign1 /= abs(sign1);
            // sign2 = lambda2 <= 0;
            // // sign2 /= abs(sign2);
            // sign3 = lambda3 <= 0;
            // // sign3 /= abs(sign3);

            // wave strength
            alpha1 = 0.5 * (pDelta - dRoe[i] * aRoe[i] * uDelta) / pow(aRoe[i], 2.0);
            alpha2 = dDelta - pDelta / pow(aRoe[i], 2.0);
            alpha3 = 0.5 * (pDelta + dRoe[i] * aRoe[i] * uDelta) / pow(aRoe[i], 2.0);

            // flux
            dFluxL = m[i];
            mFluxL = m[i] * uL + pL;
            EFluxL = d[i] * uL * hL;
            dFluxR = m[i + 1];
            mFluxR = m[i + 1] * uR + pR;
            EFluxR = d[i + 1] * uR * hR;

            dFlux[i] = 0.5 * (dFluxL + dFluxR) - 0.5 * (alpha1 * fabs(lambda1) +
                                                        alpha2 * fabs(lambda2) +
                                                        alpha3 * fabs(lambda3));
            mFlux[i] = 0.5 * (mFluxL + mFluxR) - 0.5 * (alpha1 * fabs(lambda1) * lambda1 +
                                                        alpha2 * fabs(lambda2) * lambda2 +
                                                        alpha3 * fabs(lambda3) * lambda3);
            EFlux[i] = 0.5 * (EFluxL + EFluxR) - 0.5 * (alpha1 * fabs(lambda1) * (hRoe[i] - uRoe[i] * aRoe[i]) +
                                                        alpha2 * fabs(lambda2) * (0.5 * uRoe[i] * uRoe[i]) +
                                                        alpha3 * fabs(lambda3) * (hRoe[i] + uRoe[i] * aRoe[i]));
        }

        //
    }

    /* time advancing */
    void Advancing(double dt)
    {
        for (int i = 1; i < ncell - 1; i++)
        {
            d[i] -= dt / dx * (dFlux[i] - dFlux[i - 1]);
            m[i] -= dt / dx * (mFlux[i] - mFlux[i - 1]);
            E[i] -= dt / dx * (EFlux[i] - EFlux[i - 1]);
        }
    }
};

struct FluxDifferenceSplitting
{
    /* data */
    Riemann r;

    Grid grid;

    /* functions */
    bool Init(int ncell, double xmin, double xmax, double dL, double uL, double pL, double dR, double uR, double pR);
    void Write(string fileName);

    /* time advancing */
    void AdvancingRP(double dt, double tmax);
};

bool FluxDifferenceSplitting::Init(int ncell, double xmin, double xmax, double dL, double uL, double pL, double dR, double uR, double pR)
{
    grid.InitCoordinate(ncell, xmin, xmax);
    grid.InitValue(dL, uL, pL, dR, uR, pR);

    return r.Init(dL, uL, pL, dR, uR, pR);
}

void FluxDifferenceSplitting::Write(string fileName)
{
    grid.Write(fileName);
}

void FluxDifferenceSplitting::AdvancingRP(double dt, double tmax)
{
    int step = tmax / dt;

    for (int i = 0; i < step; i++)
    {
        grid.FluxRP();
        grid.Advancing(dt);
    }
}

#endif