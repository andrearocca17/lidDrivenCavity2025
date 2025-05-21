#ifndef FIELDS_H
#define FIELDS_H
#include <vector>
#include <string>
#include "Grid.h"
#include "Solution.h"
using std::vector;

class Fields
{
    public:
        Fields();
        Fields(int&, int&);
        virtual ~Fields();

        typedef vector<Fields> vec1dField;
        typedef vector<vec1dField> vectorField;
        typedef vector<Fields::vectorField> listVectorField;

        double value;
        int NI;
        int NJ;
        int NIM; 
        int NJM;
        double X,XC,Y,YC,FXE,FYN,FXP,FYP,DXPtoE, DYPtoN;
        double Se,Sn,visc,density,volume;

        void getGridInfoPassed(Fields::vectorField&,Grid&,Solution&);
        void setVectorFieldGridFeatures(Fields::listVectorField&,Grid&,Solution&);
        void copyInternalField(Fields::vectorField&, Fields::vectorField&);
        void initalizeFields(Fields::vectorField&, double&);
        void dirichletBoundary(Fields::vectorField&,std::string&, double&);
        void oscillatingValueBoundary(Fields::vectorField&, std::string&, double&, double&,double&);
        void extaraPolateZeroGrad(Fields::vectorField&,vector<double>&,vector<double>&, std::string&);
        void shiftSolutionInTime(Fields::vectorField&, Fields::vectorField&);
        void computeCellCenterPressureGrad(Fields::vectorField&,Fields::vectorField&,Fields::vectorField&);
        
        Fields::vectorField interpolatedFieldEast(Fields::vectorField&,Grid&);
        Fields::vectorField interpolatedFieldNorth(Fields::vectorField&,Grid&);
        Fields::vectorField cellFaceGradientEast(Fields::vectorField&,Grid&);
        Fields::vectorField cellFaceGradientNorth(Fields::vectorField&,Grid&);

        void computeEastMassFluxes(Fields::vectorField&,Fields::vectorField&);
        void computeNorthMassFluxes(Fields::vectorField&,Fields::vectorField&);


        void print2dfield(Fields::vectorField&);


    protected:

    private:

};

#endif // FIELDS_H
