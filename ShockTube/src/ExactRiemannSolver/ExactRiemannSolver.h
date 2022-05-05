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

#endif