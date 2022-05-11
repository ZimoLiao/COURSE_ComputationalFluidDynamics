#include "FluxDifferenceSplitting.h"

int main()
{
    FluxDifferenceSplitting fds;

    fds.Init(100, -2.0, 2.0, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1); // Sod problem
    fds.Write("Sod.plt");

    fds.AdvancingRP(0.0001, 1.0);
    fds.Write("FDS-RP.plt");
}