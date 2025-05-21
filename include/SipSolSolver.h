#ifndef SIPSOLSOLVER_H
#define SIPSOLSOLVER_H

#include "Grid.h"
#include "Fields.h"
#include "FiniteMatrix.h"
#include "Solution.h"
#include "forAllOperations.h"
#include <string>

class SipSolSolver
{



    public:
        SipSolSolver();
        SipSolSolver(int&, int&,int&);
        virtual ~SipSolSolver();

        typedef vector<SipSolSolver> Svector1d;
        typedef vector<vector<SipSolSolver> > Svector;

        double value;
        double Residual;
        double RSM;

        Fields::vectorField solve(Fields::vectorField&,Solution&,FiniteMatrix::finiteMat& ,FiniteMatrix::finiteMat& ,FiniteMatrix::finiteMat& ,
                   FiniteMatrix::finiteMat& ,FiniteMatrix::finiteMat& ,FiniteMatrix::finiteMat&,FiniteMatrix::finiteMat&,int&,int&);

        Fields::vectorField solve(Fields::vectorField&,Solution&,FiniteMatrix::finiteMat&,FiniteMatrix::finiteMat&,int&,int&);

        void printSvector(SipSolSolver::Svector&);
    private:

            int NI, NJ, NIM, NJM, Literations;
            Svector UE,UN,LW,LS,LPR,RES;
            vector<double> RESOR;
                vector<string> EqnNames;
                vector<double> SOR;


};

#endif // SIPSOLSOLVER_H
