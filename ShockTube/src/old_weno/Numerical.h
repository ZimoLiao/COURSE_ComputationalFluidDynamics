#ifndef NUMERICAL_H_
#define NUMERICAL_H_

#include <iostream>
#include <iomanip>

using std::cin;
using std::cout;
using std::endl;

#include "Grid.h"

int factorial(int value, int order)
{
    int ans = 1;
    for (int i = 0; i != order; i++)
    {
        if (value != i)
        {
            ans *= value - i;
        }
    }
    return ans;
}

int factorial(int value)
{
    int ans = 1;
    for (int i = 0; i != value; i++)
    {
        ans *= value - i;
    }
    return ans;
}

int termial(int value, int order)
{
    int ans = 0;
    for (int i = 0; i != order; i++)
    {
        if (value != i)
        {
            ans += value - i;
        }
    }
    return ans;
}

int termial(int value)
{
    int ans = 0;
    for (int i = 0; i != value; i++)
    {
        ans += value - i;
    }
    return ans;
}

/* Newton divided differences */
class DividedDiff
{
private:
    /* data */
    int order;
    int ncell, nsize;

    double *data;

    /* functions */
    double &V(double left, double right);

public:
    DividedDiff(Grid &x, Grid &v_mean, int order);
    ~DividedDiff();

    double &operator()(double left, double right);
};

DividedDiff::DividedDiff(Grid &x, Grid &v_mean, int order)
{
#ifdef _DEBUG
    if (x.GetNcell() != v_mean.GetNcell())
    {
        cout << "divided difference | size_mismatch\n";
        data = new double[1];
        data[0] = 0.0;
        return;
    }
    if (order > x.GetNcell())
    {
        cout << "divided difference | order_out_of_range\n";
        data = new double[1];
        data[0] = 0.0;
        return;
    }
#endif
    ncell = x.GetNcell();
    this->order = order;
    nsize = termial(ncell, order);

    data = new double[nsize];
    for (int i = 1; i < ncell + 1; i++)
    {
        V(i - 0.5, i + 0.5) = v_mean(i);
    }
    // TODO: undivided differences for uniform grid
    for (int k = 2; k < order + 1; k++)
    {
        for (int i = 1; i < ncell - k + 2; i++)
        {
            V(i - 0.5, i + k - 0.5) = (V(i + 0.5, i + k - 0.5) - V(i - 0.5, i + k - 1.5)) / (x(i + k - 0.5) - x(i - 0.5));
        }
    }
}

DividedDiff::~DividedDiff()
{
    delete[] data;
}

double &DividedDiff::operator()(double left, double right)
{
    int k = right - left; // @lzmo: 浮点数0.5是精确的

#ifdef _DEBUG
    if (left <= 0 || int(left) + k > ncell || k < 1 || k > order)
    {
        cout << "divided difference | out_of_range\n";
        return data[0];
    }
#endif
    return data[termial(ncell, k - 1) + int(left)];
}

double &DividedDiff::V(double left, double right)
{
    int k = right - left; // @lzmo: 浮点数0.5是精确的

#ifdef _DEBUG
    if (left <= 0 || int(left) + k > ncell || k < 1 || k > order)
    {
        cout << "divided difference | out_of_range\n";
        return data[0];
    }
#endif
    return data[termial(ncell, k - 1) + int(left)];
}

/* Lagrangian interpolation */
void InterpL(double *x, int k, double xp, double *c)
{
    for (int i = 0; i != k; i++)
    {
        double num = 1.0, den = 1.0;
        for (int j = 0; j != k; j++)
        {
            if (j != i)
            {
                den *= (x[i] - x[j]);
                num *= (xp - x[j]);
            }
        }
        c[i] = num / den;
    }
    cout << '\n';
}

/* coefficient crj for uniform grid (for WENO) */
void Crj(int r, int k, char direction, double *c)
{
    switch (direction)
    {
    case 'r':
        for (int j = 0; j < k; j++)
        {
            c[j] = 0.0;
            for (int m = j + 1; m < k + 1; m++)
            {
                // numerator
                double num = 0.0;
                for (int l = 0; l < k + 1; l++)
                {
                    double prod = 1.0;
                    if (l != m)
                    {
                        for (int q = 0; q < k + 1; q++)
                        {
                            if (q != m && q != l)
                            {
                                prod *= (r - q + 1);
                            }
                        }
                        num += prod;
                    }
                }

                // denominator
                double den = 1.0;
                for (int l = 0; l < k + 1; l++)
                {
                    if (l != m)
                    {
                        den *= (m - l);
                    }
                }

                c[j] += num / den;
            }
        }
        break;

    case 'l':
        break;
    }
}

#endif