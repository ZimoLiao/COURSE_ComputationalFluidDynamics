#include "LinearAlgebra.h"

using namespace linear_algebra;
using namespace linear_algebra::direct;

int main()
{
    double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double B[3] = {1, 1, 1};
    double X[3];

    Info info;

    Gauss(3, A, B, X, info);
}