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
void Grid::Advance(double tDelta)
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
        switch (WENO)
        {
        case 5: // 5th-order WENO
        {
            double *vMean = new double[nCell + 2],
                   *vPosi = new double[nCell + 2],
                   *vNega = new double[nCell + 2];

            for (int j = 0; j < 3; j++)
            {
                for (int i = 0; i < nCell + 2; i++)
                {
                    vMean[i] = v[i][j];
                }
                Interpolate_WENO(nCell + 2, vMean, vPosi, vNega, 5);

                for (int i = 0; i < nCell + 2; i++)
                {
                    v3[i][j] = vNega[i];
                    v4[i][j] = vPosi[i];
                }
            }

            delete[] vMean;
            delete[] vPosi;
            delete[] vNega;
        }
        break;

        default:
            break;
        }
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

/* WENO interpolation */
void Grid::Interpolate_WENO(int n, double *vMean, double *vPosi, double *vNega, int k)
{
    switch (k)
    {
    case 1:

        break;
    case 3:

        break;
    case 5:
    {
        double q[5], alpha[5], beta[5], omega[5], alphaSum;

        for (int i = 2; i < n - 2; i++)
        {
            // negative
            for (int r = 0; r < k; r++)
            {
                q[r] = 0.0;
                for (int j = 0; j < k; j++)
                {
                    q[r] += c3rj[r][j] * vMean[i - r + j];
                }
            }

            beta[0] = 13. / 12. * pow(vMean[i] - 2. * vMean[i + 1] + vMean[i + 2], 2) + 0.25 * pow(3. * vMean[i] - 4. * vMean[i + 1] + vMean[i + 2], 2);
            beta[1] = 13. / 12. * pow(vMean[i - 1] - 2. * vMean[i] + vMean[i + 1], 2) + 0.25 * pow(vMean[i - 1] - vMean[i + 1], 2);
            beta[2] = 13. / 12. * pow(vMean[i - 2] - 2. * vMean[i - 1] + vMean[i], 2) + 0.25 * pow(vMean[i - 2] - 4. * vMean[i - 1] + 3. * vMean[i], 2);

            alphaSum = 0.0;
            for (int r = 0; r < k; r++)
            {
                alpha[r] = d3[r] / pow(1e-6 + beta[r], 2);
                alphaSum += alpha[r];
            }

            vNega[i] = 0.0;
            for (int r = 0; r < k; r++)
            {
                omega[r] = alpha[r] / alphaSum;
                vNega[i] += omega[r] * q[r];
            }

            // positive
            // TODO: 如果数组支持负指标更方便得多
            q[0] = 0.0;
            for (int j = 0; j < k; j++)
            {
                q[0] += c3rj[k - 1][k - 1 - j] * vMean[i + j];
            }
            for (int r = 1; r < k; r++)
            {
                q[r] = 0.0;
                for (int j = 0; j < k; j++)
                {
                    q[r] += c3rj[r - 1][j] * vMean[i - r + j];
                }
            }

            alphaSum = 0.0;
            for (int r = 0; r < k; r++)
            {
                alpha[r] = d3[k - 1 - r] / pow(1e-6 + beta[r], 2);
                alphaSum += alpha[r];
            }

            vPosi[i - 1] = 0.0;
            for (int r = 0; r < k; r++)
            {
                omega[r] = alpha[r] / alphaSum;
                vPosi[i - 1] += omega[r] * q[r];
            }
        }

        // TODO: 偷懒  其实应该降阶WENO
        vNega[0] = vMean[0];
        vNega[1] = vMean[1];
        vNega[n - 2] = vMean[n - 2];
        vPosi[0] = vMean[1];
        vPosi[n - 3] = vMean[n - 2];
        vPosi[n - 2] = vMean[n - 1];
    }
    break;
    case 7:

        break;
    }
}