#include <iostream>
#include "Numerical.h"
#include "Weno.h"

using namespace std;

constexpr double gamma = 1.4;

void Flux_LF(double rhon, double mn, double en,
             double rhop, double mp, double ep,
             double rhoh, double mh, double eh);

int main()
{
    int ncell = 100;
    Grid rho(ncell, 'c'), m(ncell, 'c'), e(ncell, 'c');
}

void Flux_LF(double rhon, double mn, double en,
             double rhop, double mp, double ep,
             double rhoh, double mh, double eh)
{
    // pressure
    double pn = 0.4 * (en - 0.5 * rhon * mn),
           pp = 0.4 * (ep - 0.5 * rhop * mp);

    double rhoa = 0.0,
           ma = max(fabs(2.0 * mn / rhon - 0.2 * rhon), fabs(2.0 * mp / rhop - 0.2 * rhop)),
           ea = max(fabs(1.4 * mn / rhon), fabs(1.4 * mp / rhop));

    rhoh = 0.5 * (mn + mp - rhoa * (rhop - rhon));
    mh = 0.5 * (mn * mn / rhon + pn + mp * mp / rhop + pp - ma * (mp - mn));
    eh = 0.5 * (mn * (en + pn) / rhon + mp * (ep + pp) / rhop - ea * (ep - en));
}