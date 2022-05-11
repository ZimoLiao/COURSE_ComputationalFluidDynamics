#include <iostream>
#include "math.h"
#include "Numerical.h"
#include "Eno.h"
#include "Weno.h"

using namespace std;

constexpr double gamma = 1.4;

double FluxD(double d, double m, double e);
double FluxM(double d, double m, double e);
double FluxE(double d, double m, double e);
void Flux_LF(Grid &db, Grid &mb, Grid &eb, Grid &dh, Grid &mh, Grid &eh);
void Flux_G(Grid &db, Grid &mb, Grid &eb, Grid &dh, Grid &mh, Grid &eh);

int main()
{
    /* Component-wise FV 1D shock tube */
    constexpr int ncell = 100;

    Grid x(ncell, 'm');
    Grid d(ncell, 'c'), db(ncell, 'f'),                  // density
        m(ncell, 'c'), mb(ncell, 'f'),                   // momentum
        e(ncell, 'c'), eb(ncell, 'f');                   // total energy
    Grid dh(ncell, 'b'), mh(ncell, 'b'), eh(ncell, 'b'); // flux

    // initial condition (Sod's problem)
    double DL = 1.0, ML = 0.0, PL = 1.0,
           DR = 0.125, MR = 0.0, PR = 0.1,
           EL = PL / (gamma - 1.0) + 0.5 * ML * ML / DL,
           ER = PR / (gamma - 1.0) + 0.5 * MR * MR / DR;
    x.InitUniformCoordinate(-0.5, 0.5);

    // initial condition (Lax's problem)
    // double DL = 0.445, ML = 0.698, PL = 3.528,
    //        DR = 0.5, MR = 0.0, PR = 0.1,
    //        EL = PL / (gamma - 1.0) + 0.5 * ML * ML / DL,
    //        ER = PR / (gamma - 1.0) + 0.5 * MR * MR / DR;
    // x.InitUniformCoordinate(-5, 5);

    for (int i = 1; i < ncell + 1; i++)
    {
        if (x(i) <= 0.0)
        {
            d(i) = DL;
            m(i) = ML;
            e(i) = EL;
        }
        else
        {
            d(i) = DR;
            m(i) = MR;
            e(i) = ER;
        }
        dh.l(i) = FluxD(d(i), m(i), e(i));
        mh.l(i) = FluxM(d(i), m(i), e(i));
        eh.l(i) = FluxE(d(i), m(i), e(i));
    }
    dh.r(ncell) = FluxD(d(ncell), m(ncell), e(ncell)); // TODO: 是否合理？
    mh.r(ncell) = FluxM(d(ncell), m(ncell), e(ncell));
    eh.r(ncell) = FluxE(d(ncell), m(ncell), e(ncell));

    x.Write("x.dat");
    // d.Write("d0.dat");
    // m.Write("m0.dat");
    // e.Write("e0.dat");

    // TVD-RK3
    double CFL = 0.2;
    double dt = 0.0001; // TODO: 考虑CFL条件
    double dx = (x(ncell) - x(1)) / (ncell - 1.0);

    int k = 3; // (2 * k - 1)-th order

    int step = 0;
    double T = 0.0;
    while (T < 0.2 && step < 50)
    {
        step++;

        Weno_JS(x, d, db, k);
        Weno_JS(x, m, mb, k);
        Weno_JS(x, e, eb, k);
        // Eno(x, d, db, k);
        // Eno(x, m, mb, k);
        // Eno(x, e, eb, k);

        // Flux_LF(db, mb, eb, dh, mh, eh);
        Flux_G(db, mb, eb, dh, mh, eh);

        // double umax = fabs((m / d).Max());
        // cout << umax << '\t';
        // dt = dx * CFL / (umax + 1);
        T += dt;

        for (int i = k + 1; i < ncell - k; i++)
        {
            d(i) -= dt / dx * (dh.r(i) - dh.l(i));
            m(i) -= dt / dx * (mh.r(i) - mh.l(i));
            e(i) -= dt / dx * (eh.r(i) - eh.l(i));
        }
    }

    d.Write("d.dat");
    m.Write("m.dat");
    e.Write("e.dat");
}

double FluxD(double d, double m, double e)
{
    return m;
}

double FluxM(double d, double m, double e)
{
    return m * m / d + (gamma - 1.0) * (e - 0.5 * d * m);
}

double FluxE(double d, double m, double e)
{
    return m / d * (gamma * e + 0.5 * (gamma - 1.0) * d * m);
}

void Flux_LF(Grid &db, Grid &mb, Grid &eb, Grid &dh, Grid &mh, Grid &eh)
{
    int ncell = db.GetNcell();

    for (int i = 1; i < ncell; i++) // TODO: consider boundaries
    {
        double dn = db.r(i), dp = db.l(i + 1), mn = mb.r(i), mp = mb.l(i + 1), en = eb.r(i), ep = eb.l(i + 1);

        dh.r(i) = 0.5 * (FluxD(dn, mn, en) + FluxD(dp, mp, ep));
        mh.r(i) = 0.5 * (FluxM(dn, mn, en) + FluxM(dp, mp, ep) - max(fabs(2.0 * mn / dn - (gamma - 1.0) * dn / 2.0), fabs(2.0 * mp / dp - (gamma - 1.0) * dp / 2.0)) * (mp - mn));
        eh.r(i) = 0.5 * (FluxE(dn, mn, en) + FluxE(dp, mp, ep) - max(fabs(gamma * mn / dn), fabs(gamma * mp / dp)) * (ep - en));
    }
}

void Flux_G(Grid &db, Grid &mb, Grid &eb, Grid &dh, Grid &mh, Grid &eh)
{
    int ncell = db.GetNcell();

    for (int i = 1; i < ncell; i++) // TODO: consider boundaries
    {
        double dn = db.r(i), dp = db.l(i + 1), mn = mb.r(i), mp = mb.l(i + 1), en = eb.r(i), ep = eb.l(i + 1);

        if (dn <= dp)
        {
            dh.r(i) = min(FluxD(dn, mn, en), FluxD(dp, mp, ep));
        }
        else
        {
            dh.r(i) = max(FluxD(dn, mn, en), FluxD(dp, mp, ep));
        }

        if (mn <= mp)
        {
            mh.r(i) = min(FluxM(dn, mn, en), FluxM(dp, mp, ep));
        }
        else
        {
            mh.r(i) = max(FluxM(dn, mn, en), FluxM(dp, mp, ep));
        }

        if (en <= ep)
        {
            eh.r(i) = min(FluxE(dn, mn, en), FluxE(dp, mp, ep));
        }
        else
        {
            eh.r(i) = max(FluxE(dn, mn, en), FluxE(dp, mp, ep));
        }
    }
}