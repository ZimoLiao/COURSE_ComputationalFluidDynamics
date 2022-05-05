#ifndef ERS_H_
#define ERS_H_

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

using std::cout;
using std::endl;
using std::max;
using std::ofstream;
using std::string;

/*
    Exact Riemann Solver
    lzmo, ustc, 2022-05-03
*/

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
    bool Init(double dL, double uL, double pL, double dR, double uR, double pR);

    // pressure positivity consition
    bool CheckPressurePositivity();

    // guess value of pressure in star-region (p_PV)
    double GuessPressureStar();
};

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

#endif