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
    void Advancing(double tDelta);

    // flux difference splitting
    void FluxEuler_RP(int WENO); // Roe-Pike flux

    void FluxCenter();                // v1, v2
    void VariablePlusMinus(int WENO); // v3, v4

    // flux vector splitting
};

#endif