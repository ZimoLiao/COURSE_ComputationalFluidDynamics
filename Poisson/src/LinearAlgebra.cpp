#include "LinearAlgebra.h"

namespace linear_algebra
{
    void Info::TimerOn()
    {
        start = clock();
    }

    void Info::TimerOff()
    {
        end = clock();
        time += double(end - start) / CLOCKS_PER_SEC;
    }

    namespace direct
    {
        void Gauss(int N, double *A, double *B, double *X, Info &INFO)
        {
            INFO.TimerOn();

            double buf;

            // forward elimination
            for (int i = 0; i < N; i++)
            {
                // exchange rows to ensure diagonal dominance (avoid swamp)
                int imax = i;
                for (int k = i + 1; k < N; k++)
                {
                    if (abs(A[Index(N, k, i)]) > abs(A[Index(N, i, i)]))
                    {
                        imax = k;
                    }
                }
                double Amax = A[Index(N, imax, i)];

                if (Amax != 0)
                {
                    if (imax != i)
                    {
                        for (int j = 0; j < N; j++)
                        {
                            buf = A[Index(N, i, j)];
                            A[Index(N, i, j)] = A[Index(N, imax, j)];
                            A[Index(N, imax, j)] = buf;

                            buf = B[i];
                            B[i] = B[imax];
                            B[imax] = buf;
                        }
                    }

                    // elimination
                    for (int k = i + 1; k < N; k++)
                    {
                        if (A[Index(N, k, i)] != 0)
                        {
                            buf = A[Index(N, k, i)] / Amax;

                            B[k] = B[k] - B[i] * buf;
                            for (int j = i; j < N; j++)
                            {
                                A[Index(N, k, j)] = A[Index(N, k, j)] - A[Index(N, i, j)] * buf;
                            }
                        }
                    }
                }
                else
                {
                    cout << "singular matrix\n";
                    INFO.TimerOff();
                    return;
                }
            }

            // back substitution
            X[N - 1] = B[N - 1] / A[Index(N, N - 1, N - 1)];
            for (int i = N - 2; i >= 0; i--)
            {
                buf = 0.0;

                for (int j = i + 1; j < N; j++)
                {
                    buf += A[Index(N, i, j)] * X[j];
                }

                X[i] = (B[i] - buf) / A[Index(N, i, i)];
            }

            INFO.TimerOff();
            INFO.step++;
        }

    } // namespace direct

    namespace iterative
    {

    } // namespace iterative

} // namespace linear_algebra
