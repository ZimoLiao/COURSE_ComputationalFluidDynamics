#ifndef RIEMANN_H_
#define RIEMANN_H_

#include <math.h>
#include <algorithm>

using std::max;

// Riemann problem (1D Euler equaiton)
struct Riemann
{
    /* constants */
    const double gamma = 1.4;

    const double g1 = (gamma - 1.0) / (2.0 * gamma),
                 g2 = (gamma + 1.0) / (2.0 * gamma),
                 g3 = (2.0 * gamma) / (gamma - 1.0),
                 g4 = 2.0 / (gamma - 1.0),
                 g5 = 2.0 / (gamma + 1.0),
                 g6 = (gamma - 1.0) / (gamma + 1.0),
                 g7 = (gamma - 1.0) / 2.0,
                 g8 = gamma - 1.0;

    const double tol = 1e-6;

    /* data */
    // density, velocity, pressure, sound speed
    double dL, uL, pL, cL, dR, uR, pR, cR;

    /* functions */
    // init left/right primitive variables
    bool Init(double dL, double uL, double pL, double dR, double uR, double pR)
    {
        this->dL = dL;
        this->uL = uL;
        this->pL = pL;
        this->dR = dR;
        this->uR = uR;
        this->pR = pR;

        cL = sqrt(gamma * pL / dL);
        cR = sqrt(gamma * pR / dR);

        return CheckPressurePositivity();
    }

    // pressure positivity consition
    bool CheckPressurePositivity()
    {
        return g4 * (cL + cR) >= (uR - uL);
    }

    // guess value of pressure in star-region (p_PV)
    double GuessPressureStar()
    {
        double pPV = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (dL + dR) * (cL + cR);
        return max(tol, pPV);
    }
};

#endif