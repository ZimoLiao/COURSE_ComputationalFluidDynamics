#ifndef WENO_H_
#define WENO_H_

#include <math.h>
#include "Numerical.h"

using std::fabs;

void WenoCoefficientDr() {}

// WENO interpolation for 1d uniform grid
void Weno_JS(Grid &x, Grid &v_mean, Grid &p, int k)
{
#ifdef _DEBUG
    if (x.GetType() != 'm' || v_mean.GetType() == 'b' || v_mean.GetType() == 'f' || p.GetType() == 'c')
    {
        cout << "eno | intput_type_error\n";
        return;
    }
#endif

    int ncell = x.GetNcell();

    // stencil coefficient calculation
    // for uniform grid only
    double *dr = new double[k], *dl = new double[k];
    double **cr = new double *[k];
    for (int r = 0; r < k; r++)
    {
        cr[r] = new double[2 * k - 1];
        for (int i = 0; i < 2 * k - 1; i++)
        {
            cr[r][i] = 0.0;
        }
    }
    {
        // local variables
        double *sx = new double[2 * k - 1];
        double *sc2k = new double[2 * k - 1];
        for (int i = 0; i < 2 * k - 1; i++)
        {
            sx[i] = i;
            sc2k[i] = 0.0;
        }

        for (int r = 0; r < k; r++)
        {
            Crj(r, k, 'r', &cr[r][k - 1 - r]);
        }
        Crj(k - 1, 2 * k - 1, 'r', sc2k);

        // dr calculation
        for (int r = k - 1; r >= 0; r--)
        {
            double residual = sc2k[k - 1 - r];
            for (int i = k - 1; i > r; i--)
            {
                residual -= dr[i] * cr[i][k - 1 - r];
            }
            dr[r] = residual / cr[r][k - 1 - r];
            dl[k - 1 - r] = dr[r];

#ifdef _DEBUG
            if (dr[r] < 0)
            {
                cout << "weno | stencil_error\n";
            }
#endif
        }

        // deallocation
        delete[] sc2k;
        delete[] sx;
    }

    // set-zero
    for (int i = 1; i < ncell + 1; i++)
    {
        p.l(i) = 0.0;
        p.r(i) = 0.0;
    }

    // WENO approximation
    double eps = 1e-6;
    double *omegar = new double[k], *alphar = new double[k],
           *omegal = new double[k], *alphal = new double[k],
           *beta = new double[k];
    for (int i = k; i < ncell - k + 1; i++)
    {
        // smooth indicator
        switch (k)
        {
        case 2: // 3rd-order
            beta[0] = pow(v_mean(i + 1) - v_mean(i), 2.0);
            beta[1] = pow(v_mean(i) - v_mean(i - 1), 2.0);
            break;
        case 3: // 5th-order
            beta[0] = 13.0 * pow(v_mean(i) - 2.0 * v_mean(i + 1) + v_mean(i + 2), 2.0) / 12.0 + pow(3.0 * v_mean(i) - 4.0 * v_mean(i + 1) + v_mean(i + 2), 2.0) / 4.0;
            beta[1] = 13.0 * pow(v_mean(i - 1) - 2.0 * v_mean(i) + v_mean(i + 1), 2.0) / 12.0 + pow(v_mean(i - 1) - v_mean(i + 1), 2.0) / 4.0;
            beta[2] = 13.0 * pow(v_mean(i - 2) - 2.0 * v_mean(i - 1) + v_mean(i), 2.0) / 12.0 + pow(v_mean(i - 2) - 4.0 * v_mean(i - 1) + 3.0 * v_mean(i), 2.0) / 4.0;
            break;
        case 4: // 7th-order
            beta[0] = v_mean(i) * (2107 * v_mean(i) - 9402 * v_mean(i + 1) + 7042 * v_mean(i + 2) - 1854 * v_mean(i + 3)) + v_mean(i + 1) * (11003 * v_mean(i + 1) - 17246 * v_mean(i + 2) + 4642 * v_mean(i + 3)) + v_mean(i + 2) * (7043 * v_mean(i + 2) - 3882 * v_mean(i + 3)) + 547 * pow(v_mean(i + 3), 2.0);
            beta[1] = v_mean(i - 1) * (547 * v_mean(i - 1) - 2522 * v_mean(i) + 1922 * v_mean(i + 1) - 494 * v_mean(i + 2)) + v_mean(i) * (3443 * v_mean(i) - 5966 * v_mean(i + 1) + 1602 * v_mean(i + 2)) + v_mean(i + 1) * (2843 * v_mean(i + 1) - 1642 * v_mean(i + 2)) + 267 * pow(v_mean(i + 2), 2.0);
            beta[2] = v_mean(i - 2) * (267 * v_mean(i - 2) - 1642 * v_mean(i - 1) + 1602 * v_mean(i) - 494 * v_mean(i + 1)) + v_mean(i - 1) * (2843 * v_mean(i - 1) - 5966 * v_mean(i) + 1992 * v_mean(i + 1)) + v_mean(i) * (3443 * v_mean(i) - 2522 * v_mean(i + 1)) + 547 * pow(v_mean(i + 1), 2.0);
            beta[3] = v_mean(i - 3) * (547 * v_mean(i - 3) - 3882 * v_mean(i - 2) + 4642 * v_mean(i - 1) - 1854 * v_mean(i)) + v_mean(i - 2) * (7043 * v_mean(i - 2) - 17246 * v_mean(i - 1) + 7042 * v_mean(i)) + v_mean(i - 1) * (11003 * v_mean(i - 1) - 9402 * v_mean(i)) + 2107 * pow(v_mean(i), 2.0);
            break;
        default:
            cout << "weno | order_out_of_range\n";
            return;
        }
        double asumr = 0.0, asuml = 0.0;
        for (int r = 0; r < k; r++)
        {
            alphar[r] = dr[r] / pow(eps + beta[r], 2.0);
            asumr += alphar[r];
            alphal[r] = dl[r] / pow(eps + beta[r], 2.0);
            asuml += alphal[r];
        }

        for (int r = 0; r < k; r++)
        {
            omegar[r] = alphar[r] / asumr;
            omegal[r] = alphal[r] / asuml;
        }

        // resonstruction
        double *vr = new double[k], *vl = new double[k];
        for (int r = 0; r < k; r++)
        {
            vr[r] = 0.0;
            vl[r] = 0.0;
            for (int j = 0; j < 2 * k; j++)
            {
                vr[r] += cr[r][j] * v_mean(i - k + 1 + j);
                if (r != 0)
                {
                    if (j < 2 * k - 1)
                    {
                        vl[r] += cr[r - 1][j + 1] * v_mean(i - k + 1 + j);
                    }
                }
                else
                {
                    vl[r] += cr[k - 1][2 * k - 1 - j] * v_mean(i - k + 1 + j);
                }
            }
            p.r(i) += vr[r] * omegar[r];
            p.l(i) += vl[r] * omegal[r];
        }

        delete[] vr;
        delete[] vl;
    }

    // central approximation for boundary values for non-'f' type p(x)
    if (p.GetType() != 'f')
    {
        for (int i = k - 1; i < ncell - k + 1; i++)
        {
            p.r(i) /= 2.0;
        }
    }

    // deallocation
    delete[] dr;
    delete[] dl;
    for (int r = 0; r < k; r++)
    {
        delete[] cr[r];
    }
    delete[] cr;
    delete[] omegar;
    delete[] alphar;
    delete[] omegal;
    delete[] alphal;
    delete[] beta;
}

#endif