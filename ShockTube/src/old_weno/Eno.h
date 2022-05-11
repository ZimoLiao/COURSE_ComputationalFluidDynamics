#ifndef ENO_H_
#define ENO_H_

#include <math.h>
#include "Numerical.h"

using std::fabs;

// ENO interpolation for 1d grid
void Eno(Grid &x, Grid &v_mean, Grid &p, int order)
{
#ifdef _DEBUG
    if (x.GetType() != 'm' || v_mean.GetType() == 'b' || v_mean.GetType() == 'f' || p.GetType() == 'c')
    {
        cout << "eno | intput_type_error\n";
        return;
    }
#endif

    int ncell = x.GetNcell();
    DividedDiff V(x, v_mean, order);

    /* print divided differences
        for (int k = 1; k < order + 1; k++)
        {
            for (int i = 1; i < ncell - k + 2; i++)
            {
                cout << V(i - 0.5, i + k - 0.5) << '\t';
            }
            cout << '\n';
        }
    */

    // ENO approximation
    {
        // set-zero
        for (int i = 1; i < ncell + 1; i++)
        {
            p.l(i) = 0.0;
            p.r(i) = 0.0;
        }
        if (p.GetType() == 'm')
        {
            for (int i = 1; i < ncell + 1; i++)
            {
                p(i) = v_mean(i);
            }
        }

        int *il = new int[order + 1]; // TODO: 不再保存每步的模板选取！
        il[0] = 0;
        for (int i = 1; i < ncell + 1; i++)
        {
            // stencil selection
            il[1] = i;
            for (int l = 2; l < order + 1; l++)
            {
                if (il[l - 1] != 1)
                {
                    if (il[l - 1] + l - 1 > ncell)
                    {
                        il[l] = il[l - 1] - 1;
                    }
                    else if (fabs(V(il[l - 1] - 1.5, il[l - 1] + l - 1.5)) < fabs(V(il[l - 1] - 0.5, il[l - 1] + l - 0.5))) // fabs (instead of abs) for double variables !!
                    {
                        il[l] = il[l - 1] - 1;
                    }
                    else
                    {
                        il[l] = il[l - 1];
                    }
                }
                else
                {
                    il[l] = 1;
                }
            }

            // interpolation for left/right boundaries
            for (int j = 1; j < order + 1; j++) // j-th order
            {
                double rsum = 0.0, lsum = 0.0;
                for (int m = 0; m < j; m++)
                {
                    double rprod = 1.0, lprod = 1.0;
                    for (int l = 0; l < j; l++)
                    {
                        if (l != m)
                        {
                            rprod *= (x(i + 0.5) - x(il[order] + l - 0.5));
                            lprod *= (x(i - 0.5) - x(il[order] + l - 0.5));
                        }
                    }
                    rsum += rprod;
                    lsum += lprod;
                }
                p.r(i) += V(il[order] - 0.5, il[order] + j - 0.5) * rsum;
                p.l(i) += V(il[order] - 0.5, il[order] + j - 0.5) * lsum;
            }
        }
        delete[] il;

        // central approximation for boundary values for non-'f' type p(x)
        if (p.GetType() != 'f')
        {
            for (int i = 1; i < ncell; i++)
            {
                p.r(i) /= 2.0;
            }
        }
    }
}

// v - type 'a'
void Eno(Grid &x, Grid &v, int order)
{
    int ncell = x.GetNcell();
    DividedDiff V(x, v, order);

    // ENO approximation
    // set-zero
    for (int i = 1; i < ncell + 1; i++)
    {
        v.l(i) = 0.0;
        v.r(i) = 0.0;
    }

    int *il = new int[order + 1]; // TODO: 不再保存每步的模板选取！
    il[0] = 0;
    for (int i = 1; i < ncell + 1; i++)
    {
        // stencil selection
        il[1] = i;
        for (int l = 2; l < order + 1; l++)
        {
            if (il[l - 1] != 1)
            {
                if (il[l - 1] + l - 1 > ncell)
                {
                    il[l] = il[l - 1] - 1;
                }
                else if (fabs(V(il[l - 1] - 1.5, il[l - 1] + l - 1.5)) < fabs(V(il[l - 1] - 0.5, il[l - 1] + l - 0.5))) // fabs (instead of abs) for double variables !!
                {
                    il[l] = il[l - 1] - 1;
                }
                else
                {
                    il[l] = il[l - 1];
                }
            }
            else
            {
                il[l] = 1;
            }
        }

        // interpolation for left/right boundaries
        for (int j = 1; j < order + 1; j++) // j-th order
        {
            double rsum = 0.0, lsum = 0.0;
            for (int m = 0; m < j; m++)
            {
                double rprod = 1.0, lprod = 1.0;
                for (int l = 0; l < j; l++)
                {
                    if (l != m)
                    {
                        rprod *= (x(i + 0.5) - x(il[order] + l - 0.5));
                        lprod *= (x(i - 0.5) - x(il[order] + l - 0.5));
                    }
                }
                rsum += rprod;
                lsum += lprod;
            }
            v.r(i) += V(il[order] - 0.5, il[order] + j - 0.5) * rsum;
            v.l(i) += V(il[order] - 0.5, il[order] + j - 0.5) * lsum;
        }
    }
    delete[] il;
}

#endif