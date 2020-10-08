/* Give conditions

TA  |<----- L=0.3m ----->|TB 
a1  |--------------------|a2 
    |  k,rho*C           |
*/
#include "Parameters.h"
#include "Solver.h"
#include <iostream>
#include <time.h>

/*
 *@description: uniform grid only for specific problem
 *@variables: number of nodals: n
 *@author: Fr13ndSDP
 *@date: 2020-10-04 22:11:52
*/
int initialize(int N, double L, double Tw, double TA, double A1, double K, Vec &T, Vec &position)
{
    for (int i=0; i<N; i++)
    {
        position[i] = L/(2*N)+L/N*i;
        T[i] = Tw-A1*(TA-Tw)/K*position[i];
    }
    return 0;
}

int main()
{
    clock_t start_t,end_t;
    start_t = clock();

    Parameters P;
    P.readParameters("controlDict");

    int  cnt=0;
    double Tw = 15, TWall = 0;
    Vec position(P.N,0),T(P.N,0),aW(P.N,0),aP(P.N,0),aE(P.N,0),Su(P.N,0);

    if (P.uS == 1)
    {
        cout << "unsteady calculate, solver:" << P.solver << endl;
        initialize(P.N, P.L, Tw, P.TA, P.A1, P.K, T, position);
        P.Coeff(aW,aP,aE);

        cout << "input a Wall temperature:"<<endl;
        cin >> TWall;

        while (Tw >= TWall)
        {
            P.getSource(Su, T);
            T = solve(aW, aP, aE, Su, P.solver, P.P);
            // calculaet wall temperature
            Tw = (2*P.K/(P.L/P.N)*T[0]+P.A1*P.TA)/(P.A1+2*P.K/(P.L/P.N));
            cnt++;
        }
        cout << "time = " << cnt*60 << " Tw = " << Tw << endl;
        P.outputFile(position, T);
    }
    else
    {
        cout << "steady calculate, solver:" << P.solver << endl;
        initialize(P.N, P.L, Tw, P.TA, P.A1, P.K, T, position);
        P.Coeff(aW, aP, aE);
        P.getSource(Su, T);
        T = solve(aW, aP, aE, Su, P.solver, P.P);
        // check out flux conservation
        double fluxA = (P.TA - T[0])/(1/P.A1 + P.L/P.N/(2*P.K));
        double fluxB = (-P.TB + T[P.N - 1])/(1/P.A2 + P.L/P.N/(2*P.K));
        cout << "flux in left: "<< fluxA << " flux in right: "<< fluxB << endl;

        P.outputFile(position, T);
    }
    end_t = clock();
    cout << "time consumption:" << (double)(end_t - start_t)/CLOCKS_PER_SEC << "s" << endl;
    return 0;
}
