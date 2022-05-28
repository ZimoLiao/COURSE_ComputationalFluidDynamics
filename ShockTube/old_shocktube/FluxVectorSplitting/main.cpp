#include "FluxVectorSplitting.h"

int main()
{
    FluxVectorSplitting fvs;

    fvs.Init(400, -0.5, 0.5, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1); // Sod problem
    fvs.Write("Sod.plt");

    fvs.AdvancingSW(0.0001, 0.14);
    fvs.Write("FVS-SW.plt");

    fvs.Init(400, -0.5, 0.5, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1);
    fvs.AdvancingVL(0.0001, 0.14);
    fvs.Write("FVS-VL.plt");
}