#include "FluxVectorSplitting.h"

int main()
{
    FluxVectorSplitting fvs;

    fvs.Init(200, -2.0, 2.0, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1); // Sod problem
    fvs.Write("init.plt");

    fvs.AdvancingSW(0.00001, 1.0);
    fvs.Write("SW.plt");

    fvs.Init(200, -2.0, 2.0, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1);
    fvs.AdvancingVL(0.00001, 1.0);
    fvs.Write("VL.plt");
}