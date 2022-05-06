#ifndef GRID_H_
#define GRID_H_

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>

using std::cout;
using std::endl;
using std::max;
using std::ofstream;
using std::string;

constexpr int STEPMAX = 100000;

/* information of iterative methods */
struct Info
{
    /* data */
    double time = 0.0, error = 1e-6;
    int step = 0;

    double iterErrorN2[STEPMAX], iterErrorNinf[STEPMAX];

    /* functions */
    Info()
    {
        for (int i = 0; i < STEPMAX; i++)
        {
            iterErrorN2[i] = 0.0;
            iterErrorNinf[i] = 0.0;
        }
    }

    void Reset()
    {
        time = 0.0;
        error = 1e-6;
        step = 0;
        for (int i = 0; i < STEPMAX; i++)
        {
            iterErrorN2[i] = 0.0;
            iterErrorNinf[i] = 0.0;
        }
    }

    void Print(string method)
    {
        cout << "Method: " << method << '\n'
             << "step = " << step << '\n'
             << "time = " << time << "\n\n";
    }

    void Write(string fileName)
    {
        ofstream fout;
        fout.open(fileName);

        fout << "variables = step error_N2 error_Ninf\n";

        for (int i = 1; i <= step; i++)
        {
            fout << i << '\t' << iterErrorN2[i] << '\t' << iterErrorNinf[i] << '\n';
        }

        fout.close();
    }

    void Start()
    {
        start = clock();
    }

    void End()
    {
        end = clock();
        time += double(end - start) / CLOCKS_PER_SEC;
    }

private:
    clock_t start, end;
};

/* finite difference discretized 2D-uniform-square grid (with width = 1)
    lzmo, ustc, 2022-05-06
    cell vertex, column-major
*/
struct Grid
{
    /* data */
    int nPoint, // number of points (including ghost nodes) in one direction
        nCell,  // number of cells in one direction
        nSize;  // size of data
    double h;   // spacing

    double *data;

    Grid();
    Grid(const int nCell);
    ~Grid();

protected:
    // column-major ordered
    // (i, j) in [1:(nCell-1), 1:(nCell-1)]
    inline int Index(int i, int j);

    // pointer
    int iP, jP;

    int Pointer(int i, int j);

    // boundary flags
    inline bool atW();
    inline bool atE();
    inline bool atS();
    inline bool atN();
    inline bool atBoundary();

public:
    /* functions */
    void Flush();
    void Init(const int nCell);
    void Init(Grid &obj);
    void InitPoissonRhs(const int nCell);

    void InitPoissonRes(Grid &obj, Grid &rhs);
    void InitPoissonRes(Grid &obj);

    void Write(string fileName, string varName);

    /* operators */
    // value assignment including restriction & interpolation
    void operator=(Grid &obj);

    void operator+=(Grid &obj);

    /* solvers */
    // Gauss–Seidel method
    void SolvePoissonGS(Grid &rhs, Info &info);
    void SolvePoissonGS(Grid &rhs, int stepmax, Info &info);

    // Multi-grid accelerated Gauss–Seidel method
    void SolvePoissonMGGS(Grid &rhs, int level, int n1, int n2, Info &info);
};

Grid::Grid()
{
    nCell = 2;
    nPoint = nCell + 1;
    nSize = nPoint * nPoint;

    h = 1.0 / nCell;

    data = new double[nSize];
    Flush();
}

Grid::Grid(const int nCell)
{
    this->nCell = nCell;
    nPoint = nCell + 1;
    nSize = nPoint * nPoint;

    h = 1.0 / nCell;

    data = new double[nSize];
    Flush();
}

Grid::~Grid()
{
    delete[] data;
}

inline int Grid::Index(int i, int j)
{
    // (i, j) in [1:(nCell-1), 1:(nCell-1)]
    return nPoint * j + i;
}

int Grid::Pointer(int i, int j)
{
    iP = i;
    jP = j;
    return Index(i, j);
}

inline bool Grid::atW()
{
    return iP == 1;
}

inline bool Grid::atE()
{
    return iP == nCell - 1;
}

inline bool Grid::atS()
{
    return jP == 1;
}

inline bool Grid::atN()
{
    return jP == nCell - 1;
}

inline bool Grid::atBoundary()
{
    return (atW() || atE() || atS() || atN());
}

/* functions */
void Grid::Flush()
{
    for (int ind = 0; ind < nSize; ind++)
    {
        data[ind] = 0.0;
    }
}

void Grid::Init(const int nCell)
{
    this->nCell = nCell;
    nPoint = nCell + 1;
    nSize = nPoint * nPoint;

    h = 1.0 / nCell;

    data = new double[nSize];
    Flush();
}

void Grid::Init(Grid &obj)
{
    if (nCell != obj.nCell)
    {
        nCell = obj.nCell;
        nPoint = nCell + 1;
        nSize = nPoint * nPoint;

        h = 1.0 / nCell;

        delete[] data;
        data = new double[nSize];
    }

    for (int ind = 0; ind < nSize; ind++)
    {
        data[ind] = obj.data[ind];
    }
}

void Grid::InitPoissonRhs(const int nCell)
{
    if (this->nCell != nCell)
    {
        this->nCell = nCell;
        nPoint = nCell + 1;
        nSize = nPoint * nPoint;

        h = 1.0 / nCell;

        delete[] data;
        data = new double[nSize];
    }
    Flush();

    for (int j = 1; j < nCell; j++)
    {
        for (int i = 1; i < nCell; i++)
        {
            double x = i * h, y = j * h;
            int ind = Pointer(i, j);

            // source term
            data[ind] = h * h * sin(x) * cos(y);

            // boundary conditions modification
            if (atE())
            {
                data[ind] -= y - sin(1.0) * cos(y) / 2.0;
            }
            else if (atS())
            {
                data[ind] -= -sin(x) / 2.0;
            }
            else if (atN())
            {
                data[ind] -= x - sin(x) * cos(1.0) / 2.0;
            }
        }
    }
}

void Grid::InitPoissonRes(Grid &obj, Grid &rhs)
{
    Init(rhs);

    int indP, indW, indE, indS, indN;
    for (int j = 1; j < nCell; j++)
    {
        for (int i = 1; i < nCell; i++)
        {
            indP = Index(i, j);
            indW = Index(i - 1, j);
            indE = Index(i + 1, j);
            indS = Index(i, j - 1);
            indN = Index(i, j + 1);

            data[indP] -= obj.data[indW] + obj.data[indE] + obj.data[indS] + obj.data[indN] - 4.0 * obj.data[indP];
        }
    }
}

void Grid::InitPoissonRes(Grid &obj)
{
    int indP, indW, indE, indS, indN;
    for (int j = 1; j < nCell; j++)
    {
        for (int i = 1; i < nCell; i++)
        {
            indP = Index(i, j);
            indW = Index(i - 1, j);
            indE = Index(i + 1, j);
            indS = Index(i, j - 1);
            indN = Index(i, j + 1);

            data[indP] = -(obj.data[indW] + obj.data[indE] + obj.data[indS] + obj.data[indN] - 4.0 * obj.data[indP]);
        }
    }
}

void Grid::Write(string fileName, string varName)
{
    ofstream fout;
    fout.open(fileName);

    fout << "variables = x y " + varName + "\n"
         << "zone        I = " << nCell - 1 << " J = " << nCell - 1 << " F = point\n";

    for (int j = 1; j < nCell; j++)
    {
        for (int i = 1; i < nCell; i++)
        {
            double x = i * h, y = j * h;

            int ind = Index(i, j);

            fout << x << '\t' << y << '\t' << data[ind] << '\n';
        }
    }

    fout.close();
}

/* operators */
void Grid::operator=(Grid &obj)
{
    if (this->nCell == obj.nCell)
    { // value assignment

        for (int ind = 0; ind < nSize; ind++)
        {
            data[ind] = obj.data[ind];
        }
    }
    else if (this->nCell == 2 * obj.nCell)
    { // interpolation from coarse grid to fine grid

        bool atI, atJ;
        int il, jl, ir, jr, ind;
        for (int j = 1; j < nCell; j++)
        {
            for (int i = 1; i < nCell; i++)
            {
                atI = !(i % 2);
                atJ = !(j % 2);

                il = i / 2;
                jl = j / 2;
                ir = il + 1;
                jr = jl + 1;

                ind = Index(i, j);

                if (atI)
                {
                    if (atJ)
                    {
                        data[ind] = obj.data[obj.Index(il, jl)];
                    }
                    else
                    {
                        data[ind] = 0.5 * (obj.data[obj.Index(il, jl)] + obj.data[obj.Index(il, jr)]);
                    }
                }
                else if (atJ)
                {
                    data[ind] = 0.5 * (obj.data[obj.Index(il, jl)] + obj.data[obj.Index(ir, jl)]);
                }
                else
                {
                    data[ind] = 0.25 * (obj.data[obj.Index(il, jl)] + obj.data[obj.Index(ir, jl)] + obj.data[obj.Index(il, jr)] + obj.data[obj.Index(ir, jr)]);
                }
            }
        }
    }
    else if (this->nCell == obj.nCell / 2 && !(obj.nCell % 2))
    { // direct restriction from fine grid to coarse grid

        for (int j = 1; j < nCell; j++)
        {
            for (int i = 1; i < nCell; i++)
            {
                data[Index(i, j)] = obj.data[obj.Index(i * 2, j * 2)];
            }
        }
    }
}

void Grid::operator+=(Grid &obj)
{
    if (this->nCell == obj.nCell)
    { // value assignment

        for (int ind = 0; ind < nSize; ind++)
        {
            data[ind] += obj.data[ind];
        }
    }
}

/* solvers */
void Grid::SolvePoissonGS(Grid &rhs, Info &info)
{
    // iteration
    info.Start();

    double error = 1.0, errorCrit = info.error;
    int step = 0;
    while (step < STEPMAX && error > errorCrit)
    {
        info.step++;
        step++;
        error = 0.0;

        int indP, indW, indE, indS, indN;
        for (int j = 1; j < nCell; j++)
        {
            for (int i = 1; i < nCell; i++)
            {
                indP = Index(i, j);
                indW = Index(i - 1, j);
                indE = Index(i + 1, j);
                indS = Index(i, j - 1);
                indN = Index(i, j + 1);

                double phi = (-rhs.data[indP] + data[indS] + data[indW] + data[indN] + data[indE]) / 4.0;

                error += pow(phi - data[indP], 2.0);
                // info.iterErrorNinf[info.step] = max(info.iterErrorNinf[info.step], fabs(phi - data[indP]));

                data[indP] = phi;
            }
        }
        info.iterErrorN2[info.step] = error;
        // error = info.iterErrorNinf[info.step];
    }

    info.End();
}

void Grid::SolvePoissonGS(Grid &rhs, int stepmax, Info &info)
{
    // iteration
    info.Start();

    double error = 1.0;
    int step = 0;
    while (step < stepmax)
    {
        info.step++;
        step++;
        error = 0.0;

        int indP, indW, indE, indS, indN;
        for (int j = 1; j < nCell; j++)
        {
            for (int i = 1; i < nCell; i++)
            {
                indP = Index(i, j);
                indW = Index(i - 1, j);
                indE = Index(i + 1, j);
                indS = Index(i, j - 1);
                indN = Index(i, j + 1);

                double phi = (-rhs.data[indP] + data[indS] + data[indW] + data[indN] + data[indE]) / 4.0;

                error += pow(phi - data[indP], 2.0);
                // info.iterErrorNinf[info.step] = max(info.iterErrorNinf[info.step], fabs(phi - data[indP]));

                data[indP] = phi;
            }
        }
        info.iterErrorN2[info.step] = error;
        // error = info.iterErrorNinf[info.step];
    }

    info.End();
}

void Grid::SolvePoissonMGGS(Grid &rhs, int level, int n1, int n2, Info &info)
{
    if (nCell % int(pow(2, level - 1)) || level > 4 || level < 2)
    {
        cout << "WRONG\n";
    }
    else
    {
        switch (level)
        {
        case 2:
        {
            Grid phi_f, phi_c, res_f, res_c;
            phi_f.Init(nCell);
            phi_c.Init(nCell / 2);
            res_f.Init(nCell);
            res_c.Init(nCell / 2);

            // iteration
            info.Start();

            double error = 1.0, errorCrit = info.error;
            int step = 0;
            while (step < STEPMAX && error > errorCrit)
            {
                step++;
                error = 0.0;

                // V-cicle
                SolvePoissonGS(rhs, n1, info);

                res_f.InitPoissonRes(*this, rhs);
                res_c = res_f; // restriction

                phi_c.Flush();
                phi_c.SolvePoissonGS(res_c, n1 + n2, info);

                phi_f = phi_c; // interpolation
                this->operator+=(phi_f);
                // for (int ind = 0; ind < nSize; ind++)
                // {
                //     data[ind] += phi_f.data[ind];
                // }

                SolvePoissonGS(rhs, n2, info);

                error = info.iterErrorN2[info.step];
            }

            info.End();
        }
        break;

        case 3:
        {
            Grid phiDown[3], phiUp[3], res[3];
            for (int l = 0; l < level; l++)
            {
                phiDown[l].Init(nCell / pow(2, l));
                phiUp[l].Init(nCell / pow(2, l));
                res[l].Init(nCell / pow(2, l));
            }

            // iteration
            info.Start();

            double error = 1.0, errorCrit = info.error;
            int step = 0;
            while (step < STEPMAX && error > errorCrit)
            {
                step++;
                error = 0.0;

                // V-cicle
                SolvePoissonGS(rhs, n1, info);

                res[0].InitPoissonRes(*this, rhs);
                res[1] = res[0];

                phiDown[1].Flush();
                phiDown[1].SolvePoissonGS(res[1], n1, info);

                res[1].InitPoissonRes(phiDown[1]);
                res[2] = res[1];

                phiDown[2].Flush();
                phiDown[2].SolvePoissonGS(res[2], n1 + n2, info);

                phiUp[1] = phiDown[2];
                phiDown[1] += phiUp[1];

                phiDown[1].SolvePoissonGS(res[1], n2, info);

                phiUp[0] = phiDown[1];
                this->operator+=(phiUp[0]);

                SolvePoissonGS(rhs, n2, info);

                error = info.iterErrorN2[info.step];
            }

            info.End();
        }
        break;

        default:
            break;
        }
    }
}

#endif