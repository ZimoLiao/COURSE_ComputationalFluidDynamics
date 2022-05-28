#include "FluxDifferenceSplitting.h"

int main()
{
    FluxDifferenceSplitting fds;

    fds.Init(400, -0.5, 0.5, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1); // Sod problem
    fds.Write("Sod.plt");

    fds.AdvancingRP(0.0001, 0.14);
    fds.Write("FDS-RP.plt");
}