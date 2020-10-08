#ifndef _SOLVER_
#define _SOLVER_
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
using namespace std;

typedef vector<double> Vec;
/*
 *@description: get max abs element in vector
 *@variables: 
 *@author: Fr13ndSDP
 *@date: 2020-10-06 16:07:07
*/
double maxabs(Vec v);
/*
 *@description: TDMA algorithm for solving 1D linear system
 *@variables: matrix entries: aW aP aE and source term: Su
 *@author: Fr13ndSDP
 *@date: 2020-10-04 22:01:04
*/
Vec TDMA(Vec aW, Vec aP, Vec aE, Vec Su);
/*
 *@description: Jacobi method
 *@variables: p is precision
 *@author: Fr13ndSDP
 *@date: 2020-10-06 12:39:49
*/
Vec jacobi(Vec aW, Vec aP, Vec aE, Vec Su, double p);
/*
 *@description: Guass-Seidel
 *@variables: p is precision
 *@author: Fr13ndSDP
 *@date: 2020-10-06 12:40:58
*/
Vec guassSeidel(Vec aW, Vec aP, Vec aE, Vec Su, double p);
/*
 *@description: choose solver
 *@variables: s decide the solver, p is the precision
 *@author: Fr13ndSDP
 *@date: 2020-10-07 16:40:27
*/
Vec solve(Vec aW, Vec aP, Vec aE, Vec Su, string s, double p);

#endif