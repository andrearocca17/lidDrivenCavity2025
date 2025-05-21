#ifndef GRID_H
#define GRID_H




#include <vector>
using std::vector;

class Grid
{
    public:
        Grid(int&, int&, double& );
        virtual ~Grid();

        vector<double> X;
        vector<double> XC;
        vector<double> Y;
        vector<double> YC;
        vector<double> FX;
        vector<double> FY;

        void setX(vector<double>&);
        void setXC(vector<double>&,vector<double>&);
        void setY(vector<double>&);
        void setYC(vector<double>&,vector<double>&);
        void setFX(vector<double>&,vector<double>&,vector<double>&);
        void setFY(vector<double>&,vector<double>&,vector<double>&);

        int pNI(){return NI;};
        int pNJ(){return NJ;};
        double getLength(){ return length; };
    protected:

    private:
        int N,M,NI,NJ,NIM,NJM;
        double length;
        double dx,dy;

};

#endif // GRID_H
