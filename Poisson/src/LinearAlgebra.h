#ifndef LINEARALGEBRA_H_
#define LINEARALGEBRA_H_

#include <iostream>
#include <time.h>

using std::cout;
using std::endl;

// linear algebra library (for square matrix only)
//  lzmo, ustc, 2022-05-02
//  matrices are stored in row-major form
namespace linear_algebra
{
    struct Info
    {
        /* data */
        int step = 0;
        double time = 0.0;
        double residual = 0.0;

    private:
        clock_t start, end;

    public:
        void TimerOn();
        void TimerOff();
    };

    // row-major
    inline int Index(int N, int i, int j)
    {
        return N * i + j;
    }

    // direct solver for linear equation systems
    namespace direct
    {
        void Gauss(int N, double *A, double *B, double *X, Info &INFO);

        void Thomas(int N, double *W, double *P, double *E, double *B, double *X, Info &INFO);
    } // namespace direct

    // iterative solver for linear equation systems
    namespace iterative
    {
        void PoissonJ(int N, double *RHS, double *X, Info &INFO);

        void PoissonGS(int N, double *RHS, double *X, Info &INFO);

        void PoissonSOR(int N, double *RHS, double *X, double OMEGA, Info &INFO);

        void PoissonLGSX(int N, double *RHS, double *X, Info &INFO);

        void PoissonMGGS(int N, double *RHS, double *X, Info &INFO);
    } // namespace iterative

} // namespace linear_algebra

#endif