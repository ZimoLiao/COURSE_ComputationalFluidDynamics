#include <iostream>
#include "Numerical.h"
#include "Eno.h"
#include "Weno.h"

using namespace std;

int main()
{
    constexpr int ncell = 10;
    Grid x(ncell, 'm');
    x.InitUniformCoordinate(-0.5, 0.5);
    Grid v(ncell, 'a');
    for (int i = 1; i < ncell + 1; i++)
    {
        if (x(i) <= 0)
        {
            v(i) = 1.0;
        }
    }

    Eno(x, v, 3);
    v.Print();
}