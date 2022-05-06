#include "LaxWendroff.h"

int main()
{
    Grid grid;

    grid.InitCoordinate(100, -2.0, 2.0);
    grid.InitValue(1.0, 0.0, 1.0, 0.125, 0.0, 0.1);

    grid.Write("Sod.plt");

    double tmax = 1.0, dt = 0.0001;
    int step = tmax / dt;

    for (int i = 0; i < step; i++)
    {
        grid.AdvancingLW(dt);
    }

    grid.Write("LW.plt");
}