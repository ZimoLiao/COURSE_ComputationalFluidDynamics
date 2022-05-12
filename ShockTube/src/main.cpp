#include "Grid.h"

int main()
{
    Grid shocktube;

    shocktube.InitGrid(400, 3, -2.0, 2.0);

    double vL[3] = {1.0, 0.0, 2.5},
           vR[3] = {0.125, 0.0, 0.25};
    shocktube.InitValue(0.0, vL, vR);

    for (int i = 0; i < 10000; i++)
    {
        shocktube.FluxEuler_RP(0);
        shocktube.Advance(0.0001);
    }
    shocktube.WriteAsciiEuler("out.plt");
}