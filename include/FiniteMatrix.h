#ifndef FINITEMATRIX_H
#define FINITEMATRIX_H
#include <vector>
#include <iostream>
#include "forAllOperations.h"
using namespace std;

#include "Fields.h"
using std::vector;


class FiniteMatrix
{
    public:
        FiniteMatrix();
        virtual ~FiniteMatrix();

        typedef vector<vector<FiniteMatrix> > finiteMat;

        double value;
        double aevalue;
        double awvalue;
        double asvalue;
        double anvalue;
        double svalue;

        FiniteMatrix::finiteMat interpolatedFieldEast(FiniteMatrix::finiteMat&,Grid& );
        FiniteMatrix::finiteMat interpolatedFieldNorth(FiniteMatrix::finiteMat&,Grid& );
        Fields::vectorField correctedFaceVelocityEast(Fields::vectorField&,Fields::vectorField&, Fields::vectorField&, FiniteMatrix::finiteMat&,Grid&);
        Fields::vectorField correctedFaceVelocityNorth(Fields::vectorField&,Fields::vectorField&, Fields::vectorField&, FiniteMatrix::finiteMat&,Grid&);


        void correctEastMassFluxes(Fields::vectorField&,Fields::vectorField&,FiniteMatrix::finiteMat&);
        void correctNorthMassFluxes(Fields::vectorField&,Fields::vectorField&,FiniteMatrix::finiteMat&);


        void print2dfield(FiniteMatrix::finiteMat&);
        void print2dfieldsource(FiniteMatrix::finiteMat&);

       friend FiniteMatrix::finiteMat operator+=(const FiniteMatrix::finiteMat&,  FiniteMatrix::finiteMat&);
       friend FiniteMatrix::finiteMat operator-=(const FiniteMatrix::finiteMat&,  FiniteMatrix::finiteMat&);
       friend FiniteMatrix::finiteMat operator+(const FiniteMatrix::finiteMat&, const FiniteMatrix::finiteMat&);
       friend FiniteMatrix::finiteMat operator-(const FiniteMatrix::finiteMat&, const FiniteMatrix::finiteMat&);
       friend FiniteMatrix::finiteMat operator&&(const FiniteMatrix::finiteMat&, const FiniteMatrix::finiteMat&);
       friend FiniteMatrix::finiteMat operator&&(const FiniteMatrix::finiteMat&, const Fields::vectorField&);
       friend FiniteMatrix::finiteMat operator*(const double&, const FiniteMatrix::finiteMat&);
       friend FiniteMatrix::finiteMat operator*( const FiniteMatrix::finiteMat&,const double&);
       friend FiniteMatrix::finiteMat operator/(const double&, const FiniteMatrix::finiteMat&);
       friend FiniteMatrix::finiteMat operator/(const FiniteMatrix::finiteMat&,const double&);



    protected:

    private:
};


#endif // FINITEMATRIX_H
