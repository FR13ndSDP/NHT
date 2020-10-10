#include "Parameters.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

typedef vector<double> Vec;

/*
 *@description: defualt values
 *@variables: 
 *@author: Fr13ndSDP
 *@date: 2020-10-07 16:54:18
*/
Parameters::Parameters()
{
    uS = 1;
    B_L = 3;
    B_R = 3;
    Q_L = 0;
    Q_R = 0;
    L = 0.3;
    N = 64;
    dt = 60;
    A1 = 6;
    A2 = 35;
    RHOC = 1.05e6;
    K = 0.85;
    TA = 20;
    TB = -10;
    solver = "TDMA";
    P = 1e-6;
}

Parameters::~Parameters()
{}

/*
 *@description: read parameters from control file
 *@variables: 
 *@author: Fr13ndSDP
 *@date: 2020-10-07 16:53:58
*/
int Parameters::readParameters(const char *filename)
{
    stringstream buffer;
    string line;
    string paramName;
    int paramValue_int = 0;
    double paramValue= 0;
    string paramStr;

    ifstream fin(filename);
    if (!fin.good())
    {
        string msg("parameters file not found");
        msg.append(filename);
        throw runtime_error(msg);
    }
    while (fin.good())
    {
        getline(fin, line);
        if (line[0] != '#')
        {
            buffer << line; //order is important
            buffer >> paramName;
            if (paramName.compare("unSteady") == 0)
            {
                buffer >> paramValue_int;
                uS = paramValue_int;
            }
            else if (paramName.compare("Boundary_L") == 0)
            {
                buffer >> paramValue_int;
                B_L = paramValue_int;
            }
            else if (paramName.compare("Boundary_R") == 0)
            {
                buffer >> paramValue_int;
                B_R = paramValue_int;
            }
            else if (paramName.compare("Q_L") == 0)
            {
                buffer >> paramValue;
                Q_L = paramValue;
            }
            else if (paramName.compare("Q_R") == 0)
            {
                buffer >> paramValue_int;
                Q_R = paramValue_int;
            }
            else if (paramName.compare("solver") == 0)
            {
                buffer >> paramStr;
                solver = paramStr;
            }
            else if (paramName.compare("L") == 0)
            {
                buffer >> paramValue;
                L = paramValue;
            }
            else if (paramName.compare("N") == 0)
            {
                buffer >> paramValue_int;
                N = paramValue_int;
            }
            else if (paramName.compare("dt") == 0)
            {
                buffer >> paramValue;
                dt = paramValue;
            }
            else if (paramName.compare("A1") == 0)
            {
                buffer >> paramValue;
                A1 = paramValue;
            }
            else if (paramName.compare("A2") == 0)
            {
                buffer >> paramValue;
                A2 = paramValue;
            }
            else if (paramName.compare("RHOC") == 0)
            {
                buffer >> paramValue;
                RHOC = paramValue;
            }
            else if (paramName.compare("K") == 0)
            {
                buffer >> paramValue;
                K = paramValue;
            }
            else if (paramName.compare("TA") == 0)
            {
                buffer >> paramValue;
                TA = paramValue;
            }
            else if (paramName.compare("TB") == 0)
            {
                buffer >> paramValue;
                TB = paramValue;
            }
            else if (paramName.compare("P") == 0)
            {
                buffer >> paramValue;
                P = paramValue;
            }
            else
            {
                throw runtime_error(string("unknown parameter: ").append(paramName));
            }
        }
    }
    fin.close();
    return 0;
}

/*
    *@description: output as a tecplot file
    *@variables: 
    *@author: Fr13ndSDP
    *@date: 2020-10-05 18:52:59
    */
int Parameters::outputFile(Vec vec_1, Vec vec_2)
{
    ofstream outfile;
    outfile.open("temperature.dat");
//    outfile << "Title = \"The temperature distribution\"" << endl;
//    outfile << "Variables = \"x/m\", \"T/'C\""<<endl;
    for (int i=0; i< vec_1.size(); i++)
    {
        outfile << vec_1[i] << " " << vec_2[i] << endl;
    }
    outfile.close();
    return 0;
}

/*
 *@description: Get the coefficient of internal nodals
 *@variables: 
 *@author: Fr13ndSDP
 *@date: 2020-10-04 22:42:11
*/

void Parameters::Coeff(Vec &aW, Vec &aP, Vec &aE)
{
    double dx = L/N;
    for (int i=0; i<N; i++)
    {
        aW[i] = K/dx;
        aE[i] = aW[i];
        aP[i] = aW[i]+aE[i]+ uS*RHOC*dx/dt;
    }
    aE[N-1] = 0;
    aW[0] = 0;
    // left boundary
    if (B_L == 1)
        aP[0] += K/dx;
    else if (B_L == 2)
        aP[0] += -K/dx;
    else
        aP[0] += (-K/dx + 1/(1/A1+dx/(2*K)));
    // right boundary
    if (B_R == 1)
        aP[N-1] += K/dx;
    else if (B_R == 2)
        aP[N-1] += -K/dx;
    else
        aP[N-1] += (-K/dx + 1/(1/A2+dx/(2*K)));
}

/*
 *@description: Calculate Source term
 *@variables: 
 *@author: Fr13ndSDP
 *@date: 2020-10-05 15:26:22
*/
void Parameters::getSource(Vec &Su, Vec &T)
{
    double dx = L/N;
    for (int i=1; i<N-1; i++)
    {
        Su[i] = uS*RHOC*dx/dt*T[i];
    }
    // left boundary
    if (B_L == 1)
    {

        Su[0] = uS*RHOC*dx/dt*T[0] + 2*K/dx*TA;
    }
    else if (B_L == 2)
        Su[0] = uS*RHOC*dx/dt*T[0] - Q_L;
    else
        Su[0] = uS*RHOC*dx/dt*T[0] + TA/(1/A1+dx/(2*K));
    // right boundary
    if (B_R == 1)
        Su[N-1] = uS*RHOC*dx/dt*T[N-1] + 2*K/dx*TB;
    else if (B_R == 2)
        Su[N-1] = uS*RHOC*dx/dt*T[N-1] - Q_R;
    else
        Su[N-1] = uS*RHOC*dx/dt*T[N-1] + TB/(1/A2+dx/(2*K));
}
