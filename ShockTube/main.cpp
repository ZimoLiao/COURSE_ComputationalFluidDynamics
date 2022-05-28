#include "Grid.h"

int main()
{
    Grid shocktube;

    shocktube.InitGrid(200, 3, -0.5, 0.5);

    double vL[3] = {1.0, 0.0, 2.5},
           vR[3] = {0.125, 0.0, 0.25};
    shocktube.InitValue(0.0, vL, vR);

    for (int i = 0; i < floor(0.14/0.0001); i++)
    {
        shocktube.FluxEuler_RP(5);
        shocktube.Advance(0.0001);
    }
    shocktube.WriteAsciiEuler("WENO-RP.plt");
}