#ifndef ERS_H_
#define ERS_H_

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "../Riemann.h"

using std::cout;
using std::endl;
using std::max;
using std::ofstream;
using std::string;

/*
    Exact Riemann Solver
    lzmo, ustc, 2022-05-03
*/

struct ExactRiemannSolver
{
    /* data */
    Riemann r;

    double uStar, pStar;

    /* functions */
    bool Init(double dL, double uL, double pL, double dR, double uR, double pR);

    void Solve(int n, double *s, double *dS, double *uS, double *pS);
    void Write(string fileName, int n, double *s, double *dS, double *uS, double *pS);

private:
    // pressure function
    double fL, dfL, fR, dfR;

    /* functions */
    void SolveStarRegion();
    void Sample(double s, double &dS, double &uS, double &pS);
    void UpdatePressureFunction(double p);
    void PressureFunction(double &f, double &df, double p, double dK, double pK, double cK);
};

/* ExactRiemannSolver */

bool ExactRiemannSolver::Init(double dL, double uL, double pL, double dR, double uR, double pR)
{
    return r.Init(dL, uL, pL, dR, uR, pR);
}

void ExactRiemannSolver::Solve(int n, double *s, double *dS, double *uS, double *pS)
{
    SolveStarRegion();

    for (int i = 0; i < n; i++)
    {
        Sample(s[i], dS[i], uS[i], pS[i]);
    }
}

void ExactRiemannSolver::Write(string fileName, int n, double *s, double *dS, double *uS, double *pS)
{
    ofstream fout;
    fout.open(fileName);

    fout << "variables = s d u p e\n"; // including specific internal energy

    for (int i = 0; i < n; i++)
    {
        fout << s[i] << '\t' << dS[i] << '\t' << uS[i] << '\t' << pS[i] << '\t' << pS[i] / dS[i] / (r.gamma - 1.0) << '\n';
    }

    fout.close();
}

void ExactRiemannSolver::SolveStarRegion()
{
    // initial guess
    double pOld;
    pOld = r.GuessPressureStar();
    // cout << "initial guess\n\ttpGuess " << pOld << '\n';

    // Newton–Raphson iteration
    double uDiff = r.uR - r.uL;

    double tol = r.tol, cha = 2 * tol;
    int step = 0;
    while (cha > tol && step < 10)
    {
        step++;
        // cout << "iteration step " << step << '\n';

        UpdatePressureFunction(pOld);
        pStar = max(pOld - (fL + fR + uDiff) / (dfL + dfR), tol);
        cha = 2.0 * fabs((pStar - pOld) / (pStar + pOld));

        pOld = pStar;

        // cout << "\tpStar " << pStar << "\tchange " << cha << '\n';
    }

    uStar = 0.5 * (r.uL + r.uR + fR - fL);
}

void ExactRiemannSolver::Sample(double s, double &dS, double &uS, double &pS)
{
    double dL = r.dL, uL = r.uL, pL = r.pL, cL = r.cL, dR = r.dR, uR = r.uR, pR = r.pR, cR = r.cR;

    if (s <= uStar)
    { // left of the contact discontinuity
        if (pStar <= pL)
        { // left rarefaction
            double shl = uL - cL;

            if (s <= shl)
            {
                dS = dL;
                uS = uL;
                pS = pL;
            }
            else
            {
                double stl = uStar - cL * pow((pStar / pL), r.g1);

                if (s >= stl)
                {
                    dS = dL * pow(pStar / pL, 1.0 / r.gamma);
                    uS = uStar;
                    pS = pStar;
                }
                else
                {
                    uS = r.g5 * (cL + r.g7 * uL + s);
                    double cS = r.g5 * (cL + r.g7 * (uL - s));
                    dS = dL * pow(cS / cL, r.g4);
                    pS = pL * pow(cS / cL, r.g3);
                }
            }
        }
        else
        { // left shock
            double pSL = pStar / pL;

            double sl = uL - cL * sqrt(r.g2 * pSL + r.g1);

            if (s <= sl)
            {
                dS = dL;
                uS = uL;
                pS = pL;
            }
            else
            {
                dS = dL * (pSL + r.g6) / (pSL * r.g6 + 1.0);
                uS = uStar;
                pS = pStar;
            }
        }
    }
    else
    { // right of the contact discontinuity
        if (pStar >= pR)
        { // right shock
            double pSR = pStar / pR;

            double sr = uR + cR * sqrt(r.g2 * pSR + r.g1);

            if (s >= sr)
            {
                dS = dR;
                uS = uR;
                pS = pR;
            }
            else
            {
                dS = dR * (pSR + r.g6) / (pSR * r.g6 + 1.0);
                uS = uStar;
                pS = pStar;
            }
        }
        else
        { // right rarefaction
            double shr = uR + cR;

            if (s >= shr)
            {
                dS = dR;
                uS = uR;
                pS = pR;
            }
            else
            {
                double str = uStar + cR * pow((pStar / pR), r.g1);

                if (s <= str)
                {
                    dS = dR * pow((pStar / pR), 1.0 / r.gamma);
                    uS = uStar;
                    pS = pStar;
                }
                else
                {
                    uS = r.g5 * (-cR + r.g7 * uR + s);
                    double cS = r.g5 * (-cR + r.g7 * (uR - s));
                    dS = dR * pow((cS / cR), r.g4);
                    pR = pR * pow((cS / cR), r.g3);
                }
            }
        }
    }
}

void ExactRiemannSolver::UpdatePressureFunction(double p)
{
    PressureFunction(fL, dfL, p, r.dL, r.pL, r.cL);
    PressureFunction(fR, dfR, p, r.dR, r.pR, r.cR);
}

void ExactRiemannSolver::PressureFunction(double &f, double &df, double p, double dK, double pK, double cK)
{
    if (p <= pK) // rarefaction wave
    {
        double pratio = p / pK;
        f = r.g4 * cK * (pow(pratio, r.g1) - 1.0);
        df = (1.0 / (dK * cK)) * pow(pratio, -r.g2);
    }
    else // shock wave
    {
        double aK = r.g5 / dK, bK = r.g6 * pK, q = sqrt(aK / (bK + p));
        f = (p - pK) * q;
        df = (1.0 - 0.5 * (p - pK) / (bK + p)) * q;
    }
}

#endif