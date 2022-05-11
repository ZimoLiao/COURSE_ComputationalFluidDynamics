#include "Grid.h"

/* constructor & destructor */
Grid::Grid()
{
    nCell = 1;
    nVar = 3;

    // allocation
    x = new double[nCell + 2];
    v = new double *[nCell + 2];
    flux = new double *[nCell + 2];
    v1 = new double *[nCell + 2];
    v2 = new double *[nCell + 2];
    v3 = new double *[nCell + 2];
    v4 = new double *[nCell + 2];

    for (int i = 0; i < nCell + 2; i++)
    {
        v[i] = new double[nVar];
        flux[i] = new double[nVar];
        v1[i] = new double[nVar];
        v2[i] = new double[nVar];
        v3[i] = new double[nVar];
        v4[i] = new double[nVar];
    }
}

Grid::~Grid()
{
    // deallocation
    delete[] x;
    for (int i = 0; i < nCell + 2; i++)
    {
        delete[] v[i];
        delete[] flux[i];
        delete[] v1[i];
        delete[] v2[i];
        delete[] v3[i];
        delete[] v4[i];
    }
    delete[] v;
    delete[] flux;
    delete[] v1;
    delete[] v2;
    delete[] v3;
    delete[] v4;
}

/* problem initialization */
void Grid::InitGrid(int nCell, int nVar, double xMin, double xMax)
{
    // deallocation
    delete[] x;
    for (int i = 0; i < this->nCell + 2; i++)
    {
        delete[] v[i];
        delete[] flux[i];
        delete[] v1[i];
        delete[] v2[i];
        delete[] v3[i];
        delete[] v4[i];
    }
    delete[] v;
    delete[] flux;
    delete[] v1;
    delete[] v2;
    delete[] v3;
    delete[] v4;

    this->nCell = nCell;
    this->nVar = nVar;
    this->xMin = xMin;
    this->xMax = xMax;

    // grid spacing
    xDelta = (xMax - xMin) / double(nCell);

    // allocation
    x = new double[nCell + 2];
    v = new double *[nCell + 2];
    flux = new double *[nCell + 2];
    v1 = new double *[nCell + 2];
    v2 = new double *[nCell + 2];
    v3 = new double *[nCell + 2];
    v4 = new double *[nCell + 2];

    for (int i = 0; i < nCell + 2; i++)
    {
        x[i] = xDelta * (i - 0.5) + xMin;

        v[i] = new double[nVar];
        flux[i] = new double[nVar];
        v1[i] = new double[nVar];
        v2[i] = new double[nVar];
        v3[i] = new double[nVar];
        v4[i] = new double[nVar];
    }
}

void Grid::InitValue(double xC, double *vL, double *vR)
{
    for (int i = 0; i < nCell + 2; i++)
    {
        // conservative variables
        if (x[i] < xC)
        { // left
            for (int j = 0; j < nVar; j++)
            {
                v[i][j] = vL[j];
            }
        }
        else
        { // right
            for (int j = 0; j < nVar; j++)
            {
                v[i][j] = vR[j];
            }
        }

        // other variables
        for (int j = 0; j < nVar; j++)
        {
            flux[i][j] = 0.0;
            v1[i][j] = 0.0;
            v2[i][j] = 0.0;
            v3[i][j] = 0.0;
            v4[i][j] = 0.0;
        }
    }
}

/* I/O */
void Grid::WriteAscii(string fileName)
{
    ofstream fout;
    fout.open(fileName);

    fout << "variables = x ";
    for (int j = 0; j < nVar; j++)
    {
        fout << "v" + std::to_string(j + 1) << " ";
    }
    fout << '\n';

    for (int i = 0; i < nCell + 2; i++)
    {
        fout << x[i];
        for (int j = 0; j < nVar; j++)
        {
            fout << '\t' << v[i][j];
        }
        fout << '\n';
    }

    fout.close();
}

void Grid::WriteAsciiEuler(string fileName)
{
    ofstream fout;
    fout.open(fileName);

    fout << "variables = x d u p e\n"; // including specific internal energy

    double u, du2, p, e;
    for (int i = 0; i < nCell + 2; i++)
    {
        u = v[i][1] / v[i][0];
        du2 = v[i][1] * u;
        p = (gamma - 1.0) * (v[i][2] - 0.5 * du2);
        e = p / v[i][0] / (gamma - 1.0);

        fout << x[i] << '\t' << v[i][0] << '\t' << u << '\t' << p << '\t' << e << '\n';
    }

    fout.close();
}

/* functions */
// time advancing
void Grid::Advancing(double tDelta)
{
    double tx = tDelta / xDelta;

    for (int i = 1; i < nCell + 1; i++)
    {
        for (int j = 0; j < nVar; j++)
        {
            v[i][j] -= tx * (flux[i][j] - flux[i - 1][j]);
        }
    }
}

// flux difference splitting
void Grid::FluxEuler_RP(int WENO)
{
    double dRoe, uRoe, hRoe, aRoe;

    // local variables
    double dLsr, uL, eL, pL, hL, dRsr, uR, eR, pR, hR;
    double dDelta, uDelta, pDelta;
    double alpha[3], lambda[3], coef[3];

    VariablePlusMinus(WENO);
    FluxCenter();

    for (int i = 0; i < nCell + 1; i++)
    {
        dLsr = sqrt(v3[i][0]);
        dRsr = sqrt(v4[i][0]);

        uL = v3[i][1] / v3[i][0];
        uR = v4[i][1] / v4[i][0];

        pL = (gamma - 1.0) * (v3[i][2] - 0.5 * v3[i][1] * uL);
        pR = (gamma - 1.0) * (v4[i][2] - 0.5 * v4[i][1] * uR);

        hL = (v3[i][2] + pL) / v3[i][0];
        hR = (v4[i][2] + pR) / v4[i][0];

        dDelta = v4[i][0] - v3[i][0];
        uDelta = uR - uL;
        pDelta = pR - pL;

        // Roe average values
        dRoe = dLsr * dRsr;
        uRoe = (dLsr * uL + dRsr * uR) / (dLsr + dRsr);
        hRoe = (dLsr * hL + dRsr * hR) / (dLsr + dRsr);
        aRoe = sqrt((gamma - 1.0) * (hRoe - 0.5 * uRoe * uRoe));

        // eigenvalues
        lambda[0] = uRoe - aRoe;
        lambda[1] = uRoe;
        lambda[2] = uRoe + aRoe;

        // wave strength
        alpha[0] = 0.5 * (pDelta - dRoe * aRoe * uDelta) / (aRoe * aRoe);
        alpha[1] = dDelta - pDelta / (aRoe * aRoe);
        alpha[2] = 0.5 * (pDelta + dRoe * aRoe * uDelta) / (aRoe * aRoe);

        // flux
        for (int j = 0; j < 3; j++)
        {
            coef[j] = 0.5 * alpha[j] * fabs(lambda[j]);
        }

        flux[i][0] = 0.5 * (v1[i][0] + v2[i][0]) - (coef[0] + coef[1] + coef[2]);
        flux[i][1] = 0.5 * (v1[i][1] + v2[i][1]) - (coef[0] * lambda[0] + coef[1] * lambda[1] + coef[2] * lambda[2]);
        flux[i][2] = 0.5 * (v1[i][2] + v2[i][2]) - (coef[0] * (hRoe - uRoe * aRoe) + coef[1] * (0.5 * uRoe * uRoe) + coef[2] * (hRoe - uRoe * aRoe));
    }
}

void Grid::FluxCenter() // v1
{
    double u, du2, p;
    for (int i = 0; i < nCell + 2; i++)
    {
        // left
        u = v3[i][1] / v3[i][0];
        du2 = v3[i][1] * u;
        p = (gamma - 1.0) * (v3[i][2] - 0.5 * du2);

        v1[i][0] = v3[i][1];
        v1[i][1] = du2 + p;
        v1[i][2] = u * (v3[i][2] + p);

        // right
        u = v4[i][1] / v4[i][0];
        du2 = v4[i][1] * u;
        p = (gamma - 1.0) * (v4[i][2] - 0.5 * du2);

        v2[i][0] = v4[i][1];
        v2[i][1] = du2 + p;
        v2[i][2] = u * (v4[i][2] + p);
    }
}

void Grid::VariablePlusMinus(int WENO) // v3, v4
{
    if (WENO)
    {
    }
    else
    {
        for (int i = 0; i < nCell + 1; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                v3[i][j] = v[i][j];
                v4[i][j] = v[i + 1][j];
            }
        }
    }
}