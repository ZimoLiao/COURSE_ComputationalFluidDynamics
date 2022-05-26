#ifndef GRID_H_
#define GRID_H_

#include <string>
#include <fstream>
#include "math.h"

using std::ofstream;
using std::string;

/*
    Uniform grid for 1D hyperbolic conservation laws
    lzmo, ustc, 2022-05-11
*/
class Grid
{
private:
    /* constants */
    const double gamma = 1.4;

    /* parameters */
    // size
    int nCell, nVar;

    // geometry
    double xMin, xMax, xDelta;

    /* variables */
    double *x;                     // coordinates
    double **v;                    // conservative variables
    double **flux;                 // numerical intercell flux
    double **v1, **v2, **v3, **v4; // intermediate variables

public:
    /* constructor & destructor */
    Grid();
    ~Grid();

    /* problem initialization */
    void InitGrid(int nCell, int nVar, double xMin, double xMax);
    void InitValue(double xC, double *vL, double *vR);

    /* input & output */
    void WriteAscii(string fileName);
    void WriteAsciiEuler(string fileName);

    /* solver */

    // private:
    /* functions */
    // time advancing
    void Advance(double tDelta);

    // flux difference splitting
    void FluxEuler_RP(int WENO); // Roe-Pike flux

    void FluxCenter();                // v1, v2
    void VariablePlusMinus(int WENO); // v3, v4

    // flux vector splitting

private:
    /* WENO interpolation */
    void Interpolate_WENO(int n, double *vMean, double *vPosi, double *vNega, int k);

    // interpolation constants ckrj
    double c1rj = 1,
           c2rj[2][2] = {0.5, 0.5,
                         -0.5, 1.5},
           c3rj[3][3] = {1. / 3., 5. / 6., -1. / 6.,
                         -1. / 6., 5. / 6., 1. / 3.,
                         1. / 3., -7. / 6., 11. / 6.},
           c4rj[4][4] = {0.25, 13. / 12., -5. / 12., 1. / 12.,
                         -1. / 12., 7. / 12., 7. / 12., -1. / 12.,
                         1. / 12., -5. / 12., 13. / 12., 0.25,
                         -0.25, 13. / 12., -23. / 12., 25. / 12.};

    double d1 = 1,
           d2[2] = {2. / 3., 1. / 3.},
           d3[3] = {0.3, 0.6, 0.1},
           d4[4] = {1. / 35., 12. / 35., 18. / 35., 4. / 35.};
};

#endif