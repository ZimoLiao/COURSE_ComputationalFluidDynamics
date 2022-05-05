#ifndef DATA_H_
#define DATA_H_

#include <iostream>
#include <string>
#include <fstream>
#include <math.h>

using std::cout;
using std::ofstream;
using std::string;

// index calculation (colume-major)
inline int Index(int i, int j, int num)
{
    return (num + 1) * j + i;
}

// finite difference discretized 2D-uniform-square grid (x,y) in [0,1;0,1]
// cell-vertex
struct Grid
{
    /* data */
    //  num     number of cells
    //  size    size of data
    int num = 1, size = 1;

    double *data;

    Grid()
    {
        data = new double[1];
        data[0] = 0.0;
    }

    Grid(int num)
    {
        this->num = num;

        size = pow(num + 1, 2);
        data = new double[size];
        for (int ind = 0; ind < size; ind++)
        {
            data[ind] = 0.0;
        }
    }

    ~Grid()
    {
        delete[] data;
    }

    void operator+=(const Grid &obj)
    {
        if (obj.num != num)
        {
            cout << "WRONG\n";
            return;
        }

        for (int ind = 0; ind < size; ind++)
        {
            data[ind] += obj.data[ind];
        }
    }

    void Init(int num)
    {
        this->num = num;

        size = pow(num + 1, 2);
        data = new double[size];
        for (int ind = 0; ind < size; ind++)
        {
            data[ind] = 0.0;
        }
    }

    void WriteAscii(string fname)
    {
        double h = 1.0 / num;

        ofstream fout;
        fout.open(fname);

        fout << "variables = x y phi\n"
             << "zone        I = " << num - 1 << " J = " << num - 1 << " F = point\n";

        int ind;
        double x, y;
        for (int j = 1; j < num; j++)
        {
            for (int i = 1; i < num; i++)
            {
                x = i * h;
                y = j * h;

                ind = Index(i, j, num);

                fout << x << '\t' << y << '\t' << data[ind] << '\n';
            }
        }

        fout.close();
    }

    // direct restriction
    void MGrestriction(const Grid &fine)
    {
        if (fine.num != num * 2)
        {
            cout << "WRONG\n";
            return;
        }

        int numf = fine.num;

        for (int j = 1; j < num; j++)
        {
            for (int i = 1; i < num; i++)
            {
                data[Index(i, j, num)] = fine.data[Index(i * 2, j * 2, numf)];
            }
        }
    }

    // bi-linear interpolaiton
    // boundary value = zero
    void MGinterpolation(const Grid &coarse)
    {
        if (coarse.num != num / 2)
        {
            cout << "WRONG\n";
            return;
        }

        int numc = coarse.num;

        for (int j = 1; j < num; j++)
        {
            for (int i = 1; i < num; i++)
            {
                bool ati = !(i % 2), atj = !(j % 2);

                int il = i / 2, jl = j / 2,
                    ir = il + 1, jr = jl + 1;

                int index = Index(i, j, num);

                if (ati)
                {
                    if (atj)
                    {
                        data[index] = coarse.data[Index(il, jl, numc)];
                    }
                    else
                    {
                        data[index] = 0.5 * (coarse.data[Index(il, jl, numc)] + coarse.data[Index(il, jl + 1, numc)]);
                    }
                }
                else if (atj)
                {
                    data[index] = 0.5 * (coarse.data[Index(il, jl, numc)] + coarse.data[Index(il + 1, jl, numc)]);
                }
                else
                {
                    data[index] = 0.25 * (coarse.data[Index(il, jl, numc)] + coarse.data[Index(il, jl + 1, numc)] + coarse.data[Index(il + 1, jl, numc)] + coarse.data[Index(il + 1, jl + 1, numc)]);
                }
            }
        }
    }
};

struct Rhs
{
    /* data */
    int num, size;

    double *data;

    // default r.h.s. for Poisson equation (including given B.C.)
    Rhs(int num)
    {
        this->num = num;
        size = pow(num + 1, 2);
        data = new double[size];

        double h = 1.0 / num;

        int ind;
        double x, y;
        for (int j = 1; j < num; j++)
        {
            for (int i = 1; i < num; i++)
            {
                x = i * h;
                y = j * h;

                ind = Index(i, j, num);

                // source term
                data[ind] = h * h * sin(x) * cos(y);

                // boundary conditions
                if (i == num - 1)
                {
                    data[ind] -= y - sin(1.0) * cos(y) / 2.0;
                }
                if (j == 1)
                {
                    data[ind] -= -sin(x) / 2.0;
                }
                else if (j == num - 1)
                {
                    data[ind] -= x - sin(x) * cos(1.0) / 2.0;
                }
            }
        }
    }

    Rhs(const Rhs &obj)
    {
        num = obj.num;
        size = pow(num + 1, 2);
        data = new double[size];

        for (int ind = 0; ind < size; ind++)
        {
            data[ind] = obj.data[ind];
        }
    }

    ~Rhs()
    {
        delete[] data;
    }

    void Init(int num)
    {
        this->num = num;
        size = pow(num + 1, 2);
        data = new double[size];

        double h = 1.0 / num;

        int ind;
        double x, y;
        for (int j = 1; j < num; j++)
        {
            for (int i = 1; i < num; i++)
            {
                x = i * h;
                y = j * h;

                ind = Index(i, j, num);

                // source term
                data[ind] = h * h * sin(x) * cos(y);

                // boundary conditions
                if (i == num - 1)
                {
                    data[ind] -= y - sin(1.0) * cos(y) / 2.0;
                }
                if (j == 1)
                {
                    data[ind] -= -sin(x) / 2.0;
                }
                else if (j == num - 1)
                {
                    data[ind] -= x - sin(x) * cos(1.0) / 2.0;
                }
            }
        }
    }

    void Residual(Grid &grid)
    {
        for (int j = 1; j < num; j++)
        {
            for (int i = 1; i < num; i++)
            {
                int indp = Index(i, j, num);
                int inds = Index(i, j - 1, num);
                int indw = Index(i - 1, j, num);
                int indn = Index(i, j + 1, num);
                int inde = Index(i + 1, j, num);

                data[indp] -= grid.data[inds] + grid.data[indw] + grid.data[indn] + grid.data[inde] - 4.0 * grid.data[indp];
            }
        }
    }

    // direct restriction
    void MGrestriction(const Rhs &fine)
    {
        if (fine.num != num * 2)
        {
            cout << "WRONG\n";
            return;
        }

        int numf = fine.num;

        for (int j = 1; j < num; j++)
        {
            for (int i = 1; i < num; i++)
            {
                data[Index(i, j, num)] = fine.data[Index(i * 2, j * 2, numf)];
            }
        }
    }
};

#endif