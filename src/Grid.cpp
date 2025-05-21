#include "Grid.h"
#include "forAllOperations.h"
Grid::Grid(int& N_, int& M_, double& length_)
: N(N_),M(M_),length(length_), X(N_+2,0.0),XC(N_+2,0.0),Y(M_+2,0.0),YC(M_+2,0.0)
,FX(N_+2,0.0),FY(M_+2,0.0)
{
    NI=N+2; // MAX SIZE
    NJ=M+2;
    NIM=NI-1;
    NJM=NJ-1;
    dx=length/N;
    dy=length/M;
    //ctor

    setX(X);
    setXC(X,XC);
    setY(Y);
    setYC(Y,YC);
    setFX(FX,X,XC);
    setFY(FY,Y,YC);

}

Grid::~Grid()
{
    //dtor
}

void Grid::setX(vector<double>& vecX)
{

        for(int i=1; i<vecX.size(); i++)
        {
            vecX[i] =vecX[i-1] + dx;
         //   cout << vecX[i] << ' ';
        }

         vecX[vecX.size()-1] =     vecX[vecX.size()-2];
       // vecX[NI-1] = vecX[NI-2];

}


void Grid::setXC(vector<double>& vecX,vector<double>& vecXC)
{

 for(int i=1; i<vecXC.size(); i++)
        {
            vecXC[i] =0.5* (vecX[i]+vecX[i-1]);
         //   cout << vecX[i] << ' ';
        }
        vecXC[0] = vecX[0];
        vecXC[vecX.size()-1] = vecX[vecX.size()-2];


}

void Grid::setY(vector<double>& vecY)
{
        for(int i=1; i<vecY.size(); i++)
        {
            vecY[i] =vecY[i-1] + dy;
         //   cout << vecX[i] << ' ';
        }

         vecY[vecY.size()-1] =     vecY[vecY.size()-2];

}
void Grid::setYC(vector<double>& vecY,vector<double>& vecYC)
{
 for(int i=1; i<vecY.size(); i++)
        {
            vecYC[i] =0.5* (vecY[i]+vecY[i-1]);
         //   cout << vecX[i] << ' ';
        }
        vecYC[0] = vecY[0];
        vecYC[vecYC.size()-1] = vecY[vecY.size()-2];


}


void Grid::setFX(vector<double>& result,vector<double>& vecX,vector<double>& vecXC)
{
    //ERRORE QUI. 
        for(int i=0; i<vecX.size()-1; i++)
        {
            result[i] =(vecX[i]-vecXC[i])/(vecXC[i+1]-vecXC[i]);
         //   cout << vecX[i] << ' ';
        }
}
void Grid::setFY(vector<double>& result,vector<double>& vecY,vector<double>& vecYC)
{

        for(int i=0; i<vecY.size()-1; i++)
        {
            result[i] =(vecY[i]-vecYC[i])/(vecYC[i+1]-vecYC[i]);
         //   cout << vecX[i] << ' ';
        }

}

