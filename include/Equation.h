#ifndef EQUATION_H
#define EQUATION_H
#include <vector>
#include "Fields.h"
#include "FiniteMatrix.h"

class Equation
{
    public:
        Equation(const FiniteMatrix::finiteMat&);
        virtual ~Equation();
        typedef   vector<FiniteMatrix> Svector1d;
        typedef   FiniteMatrix::finiteMat Svector;


        void relax(Fields::vectorField&);
        void resetEqn();
        void noWallShearXBoundaryConditions(Fields::vectorField&);
        void noWallShearYBoundaryConditions(Fields::vectorField&);
        void assembleEquation();
        Fields::vectorField solve(Fields::vectorField&,FiniteMatrix::finiteMat&,Solution&,int&);


            FiniteMatrix::finiteMat AP2;

        double value;
        double Residual;
        double RSM;
        double RESOR;


        double URF;
        double rURF;
        string EqnName;
        double SOR;

        FiniteMatrix::finiteMat APinitial;
        FiniteMatrix::finiteMat AP;
        FiniteMatrix::finiteMat AE;
        FiniteMatrix::finiteMat AW;
        FiniteMatrix::finiteMat AS;
        FiniteMatrix::finiteMat AN;
        FiniteMatrix::finiteMat rAP;
        FiniteMatrix::finiteMat APU;
        FiniteMatrix::finiteMat sourceInitial;
        FiniteMatrix::finiteMat sourceB;
        FiniteMatrix::finiteMat sourceFinal;
        FiniteMatrix::finiteMat sourceRelaxed;

    protected:

    private:
        int NI, NJ, NIM, NJM, Literations;
        Svector UE,UN,LW,LS,LPR,RES;


};

#endif // EQUATION_H
