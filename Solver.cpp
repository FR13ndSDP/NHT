#include "Solver.h"

double maxabs(Vec v)
{
    double min = *min_element(v.begin(),v.end());
    double max = *max_element(v.begin(),v.end()); 
    if (abs(max) > abs(min))
        return abs(max);
    else
        return abs(min);
}

Vec TDMA(Vec aW, Vec aP, Vec aE, Vec Su)
{
    int size = Su.size();
    Vec P(size, 0), Q(size, 0), phi(size,0);
    P[0] = aE[0] / aP[0];
    Q[0] = Su[0] / aP[0];
    for (int i = 1; i < size; i++)
    {
        P[i] = aE[i]/(aP[i] - aW[i]*P[i-1]);
        Q[i] = (Su[i]+aW[i]*Q[i-1])/(aP[i] - aW[i]*P[i-1]);
    }
    phi[size-1] = Q[size-1];
    for (int i = size-1; i>0 ; i--)
    {
        phi[i-1] = P[i-1]*phi[i]+Q[i-1];
    }
    return phi;
}

Vec jacobi(Vec aW, Vec aP, Vec aE, Vec Su, double p)
{
    int n = aW.size();
    Vec phi(n,1);
    Vec phi0(n,0);
    Vec errList(n,1);
    // whith origin value as 0
    while (maxabs(errList) > p)
    {
        phi[0] = (aE[0]*phi0[1]+Su[0])/aP[0];
        for (int i =1; i<n-1; i++)
        {
            phi[i] = (aW[i]*phi0[i-1]+aE[i]*phi0[i+1]+Su[i])/aP[i];
        }
        phi[n-1] = (aW[n-1]*phi0[n-2]+Su[n-1])/aP[n-1];
        for (int i=0; i<n; i++)
        {
            errList[i] = phi0[i] - phi[i];
            phi0[i] = phi[i];
        }
    }
    return phi;
}

Vec guassSeidel(Vec aW, Vec aP, Vec aE, Vec Su, double p)
{
    int n = aW.size();
    Vec phi(n,1);
    Vec phi0(n,0);
    Vec errList(n,1);
    // whith origin value as 0
    while (maxabs(errList) > p)
    {
        phi[0] = (aE[0]*phi0[1]+Su[0])/aP[0];
        for (int i =1; i<n-1; i++)
        {
            phi[i] = (aW[i]*phi[i-1]+aE[i]*phi0[i+1]+Su[i])/aP[i];
        }
        phi[n-1] = (aW[n-1]*phi[n-2]+Su[n-1])/aP[n-1];
        for (int i =0; i<n; i++)
        {
            errList[i] = phi0[i] - phi[i];
            phi0[i] = phi[i];
        }
    }
    return phi;
}

Vec solve(Vec aW, Vec aP, Vec aE, Vec Su, string s, double p)
{
    if (s.compare("TDMA") == 0)
    {
        return TDMA(aW, aP, aE, Su);
    }
    else if (s.compare("jacobi") == 0)
    {
        return jacobi(aW, aP, aE, Su, p);
    }
    else
    {
        return guassSeidel(aW, aP, aE, Su, p);
    }
}
