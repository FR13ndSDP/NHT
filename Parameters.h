#ifndef _PARAMETER_
#define _PARAMETER_

#include <vector>
#include <string>
using namespace std;

class Parameters
{
private:
    
public:
    int  N, uS, B_L, B_R;
    double dt, L, A1, A2, RHOC, K, TA, TB, P, Q_R, Q_L;
    string solver;
    Parameters();
    ~Parameters();
    int readParameters(const char *filename);
    void Coeff(vector<double> &aW, vector<double> &aP, vector<double> &aE);
    void getSource(vector<double> &Su, vector<double> &T);
    int outputFile(vector<double> vec_1, vector<double> vec_2);
};

#endif
