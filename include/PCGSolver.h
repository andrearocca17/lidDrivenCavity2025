#ifndef PCGSOLVER_H
#define PCGSOLVER_H
#include "SipSolSolver.h"

#include <string>

class PCGSolver
{
    public:
        PCGSolver();
          PCGSolver(int&, int&, int&);
        virtual ~PCGSolver();

        typedef vector<SipSolSolver> Svector1d;
        typedef vector<vector<SipSolSolver> > Svector;

        double value;
        double Residual;
        double RSM;
        double RES0;
        int InnerIter;
        int InnerItermax;


        Fields::vectorField solve(Fields::vectorField&,Solution&,FiniteMatrix::finiteMat& ,FiniteMatrix::finiteMat& ,FiniteMatrix::finiteMat& ,
                   FiniteMatrix::finiteMat& ,FiniteMatrix::finiteMat& ,FiniteMatrix::finiteMat&,FiniteMatrix::finiteMat&,int&,int&);

        void printSvector(SipSolSolver::Svector&);

    protected:

    private:


            int NI, NJ, NIM, NJM, Literations;
            Svector D,ZK,PK,RES;
            vector<double> RESOR;
                vector<string> EqnNames;
                vector<double> SOR;
};

#endif // PCGSOLVER_H
