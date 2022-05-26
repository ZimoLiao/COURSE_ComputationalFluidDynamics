#ifndef GRID_H_
#define GRID_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using std::cin;
using std::cout;
using std::endl;
using std::string;

class Grid
{
private:
    /* data */
    //  c   cell center
    //  b   cell boundary
    //  m   cell center & boundary
    //  f   two different boundaries (for left/right flux)
    //  a   left/right boundaries & cell center
    char type;
    int ncell;

    double *center, *bound;

public:
    Grid(const Grid &grid);
    Grid(int ncell);
    Grid(int ncell, char type);
    ~Grid();

    double &operator()(int i);
    double &operator()(double i);
    double &l(int i);
    double &r(int i);
    double &plus(double i);
    double &minus(double i);

    Grid operator/(Grid &grid);

    int GetNcell();
    char GetType();
    double Max();

    void InitUniformCoordinate(double xl, double xr);

    void Print();
    void Write(const char *fname);
};

Grid::Grid(const Grid &grid)
{
    this->ncell = grid.ncell;
    this->type = grid.type;

    switch (type)
    {
    case 'c':
        center = new double[ncell];
        bound = new double[1];
        bound[0] = 0.0;

        for (int i = 0; i != ncell; i++)
            center[i] = grid.center[i];
        break;
    case 'b':
        bound = new double[ncell + 1];
        center = new double[1];
        center[0] = 0.0;

        for (int i = 0; i != ncell + 1; i++)
            bound[i] = grid.bound[i];
        break;
    case 'f':
        bound = new double[2 * ncell];
        center = new double[1];
        center[0] = 0.0;

        for (int i = 0; i != 2 * ncell; i++)
            bound[i] = grid.bound[i];
        break;
    case 'm':
        center = new double[ncell];
        bound = new double[ncell + 1];

        for (int i = 0; i != ncell; i++)
        {
            center[i] = grid.center[i];
            bound[i] = grid.bound[i];
        }
        bound[ncell] = grid.bound[ncell];
        break;
    }
}

Grid::Grid(int ncell)
{
    this->type = 'm';
    this->ncell = ncell;
    center = new double[ncell];
    bound = new double[ncell + 1];

    for (int i = 0; i != ncell; i++)
    {
        center[i] = 0.0;
        bound[i] = 0.0;
    }
    bound[ncell] = 0.0;
}

Grid::Grid(int ncell, char type)
{
    this->ncell = ncell;
    this->type = type;

    switch (type)
    {
    case 'c':
        center = new double[ncell];
        bound = new double[1];
        bound[0] = 0.0;

        for (int i = 0; i != ncell; i++)
            center[i] = 0.0;
        break;
    case 'b':
        bound = new double[ncell + 1];
        center = new double[1];
        center[0] = 0.0;

        for (int i = 0; i != ncell + 1; i++)
            bound[i] = 0.0;
        break;
    case 'f':
        bound = new double[2 * ncell];
        center = new double[1];
        center[0] = 0.0;

        for (int i = 0; i != 2 * ncell; i++)
            bound[i] = 0.0;
        break;
    case 'a':
        bound = new double[2 * ncell];
        center = new double[ncell];

        for (int i = 0; i != ncell; i++)
            center[i] = 0.0;
        for (int i = 0; i != 2 * ncell; i++)
            bound[i] = 0.0;
        break;
    default:
        this->type = 'm';
        center = new double[ncell];
        bound = new double[ncell + 1];

        for (int i = 0; i != ncell; i++)
        {
            center[i] = 0.0;
            bound[i] = 0.0;
        }
        bound[ncell] = 0.0;
        break;
    }
}

Grid::~Grid()
{
    delete[] center;
    delete[] bound;
}

double &Grid::operator()(int i)
{
#ifdef _DEBUG
    if (type == 'b')
    {
        cout << "grid | type_mismatch\n";
        return center[0];
    }

    if (i < 1 || i > ncell)
    {
        cout << "grid | center_out_of_range\n";
        return center[0];
    }
#endif

    return center[i - 1];
}

double &Grid::operator()(double i)
{
#ifdef _DEBUG
    if (type == 'c' || type == 'f' || type == 'a')
    {
        cout << "grid | type_mismatch\n";
        return bound[0];
    }

    if (i <= 0.0 || i >= ncell + 1.0)
    {
        cout << "grid | bound_out_of_range\n";
        return bound[0];
    }
#endif

    return bound[int(i)];
}

double &Grid::l(int i)
{
#ifdef _DEBUG
    if (type == 'c')
    {
        cout << "grid | type_mismatch\n";
        return center[0];
    }

    if (i < 1 || i > ncell)
    {
        cout << "grid | center_out_of_range\n";
        return bound[0];
    }
#endif

    switch (type)
    {
    case 'f':
        return bound[2 * (i - 1)];
    case 'a':
        return bound[2 * (i - 1)];
    default:
        return bound[i - 1];
    }
}

double &Grid::r(int i)
{
#ifdef _DEBUG
    if (type == 'c')
    {
        cout << "grid | type_mismatch\n";
        return center[0];
    }

    if (i < 1 || i > ncell)
    {
        cout << "grid | center_out_of_range\n";
        return bound[0];
    }
#endif

    switch (type)
    {
    case 'f':
        return bound[2 * i - 1];
    case 'a':
        return bound[2 * i - 1];
    default:
        return bound[i];
    }
}

// TODO: 整个解决方案需要优化！
double &Grid::plus(double i)
{
    return l(int(i) + 1);
}

double &Grid::minus(double i)
{
    return r(int(i));
}

// TODO: 随手加的
Grid Grid::operator/(Grid &grid)
{
    Grid ans(ncell, 'c');
    for (int i = 1; i < ncell + 1; i++)
    {
        ans(i) = center[i - 1] / grid(i);
    }
    return ans;
}

int Grid::GetNcell()
{
    return ncell;
}

char Grid::GetType()
{
    return type;
}

double Grid::Max()
{
    double data_max = center[0];
    for (int i = 0; i < ncell; i++)
    {
        if (data_max < center[i])
        {
            data_max = center[i];
        }
    }
    return data_max;
}

void Grid::InitUniformCoordinate(double xl, double xr)
{
#ifdef _DEBUG
    if (xr <= xl)
    {
        cout << "grid | range_error\n";
        return;
    }
#endif

    double h = (xr - xl) / ncell;

    switch (type)
    {
    case 'c':
        for (int i = 0; i != ncell; i++)
        {
            center[i] = xl + (0.5 + i) * h;
        }
        break;
    case 'b':
        for (int i = 0; i != ncell + 1; i++)
        {
            bound[i] = xl + i * h;
        }
        break;
    default:
        for (int i = 0; i != ncell; i++)
        {
            center[i] = xl + (0.5 + i) * h;
        }
        for (int i = 0; i != ncell + 1; i++)
        {
            bound[i] = xl + i * h;
        }
        break;
    }
}

void Grid::Print()
{
    bool pc = true, pb = true;
    switch (type)
    {
    case 'c':
        pb = false;
        break;
    case 'b':
        pc = false;
        break;
    case 'f':
        pb = false;
        pc = false;
        break;
    case 'a':
        pb = false;
        break;
    }

    if (pc)
    {
        cout << "    cell center v(i):\n"
             << std::setw(8) << ' ';
        for (int i = 0; i != ncell; i++)
        {
            cout << std::setprecision(4) << std::setw(8) << center[i] << std::setprecision(4) << std::setw(8) << ' ';
        }
        cout << '\n';
    }
    if (pb)
    {
        cout << "    cell boundary v(i+-0.5):\n";
        for (int i = 0; i != ncell + 1; i++)
        {
            cout << std::setprecision(4) << std::setw(8) << bound[i] << std::setprecision(4) << std::setw(8) << ' ';
        }
        cout << '\n';
    }
    if ((!pc && !pb) || type == 'a')
    {
        cout << "    left boundary v+ (i-0.5):\n"
             << std::setw(1) << ' ' << std::setw(1) << ' ';
        for (int i = 1; i < ncell + 1; i++)
        {
            cout << std::setprecision(4) << std::setw(8) << this->l(i) << std::setprecision(4) << std::setw(8) << ' ';
        }
        cout << "\n    right boundary v- (i+0.5):\n"
             << std::setprecision(4) << std::setw(7) << ' ' << std::setprecision(4) << std::setw(7) << ' ';
        for (int i = 1; i < ncell + 1; i++)
        {
            cout << std::setw(8) << this->r(i) << std::setw(8) << ' ';
        }
        cout << '\n';
    }
}

void Grid::Write(const char *fname)
{
    std::ofstream fout;
    fout.open(fname);

    for (int i = 0; i < ncell; i++)
    {
        fout << center[i] << '\n';
    }

    fout.close();
}

#endif