#include "PCGSolver.h"
#include <cmath>
PCGSolver::PCGSolver()
{
    //ctor
}

PCGSolver::PCGSolver(int& NI_, int& NJ_,int& Liter)
: NI(NI_),NJ(NJ_),NIM(NI-1),NJM(NJ-1),
D(NI,Svector1d(NJ)),ZK(NI,Svector1d(NJ)),PK(NI,Svector1d(NJ)),
RES(NI,Svector1d(NJ)) , value(0.0),Residual(0.0),RSM(0.0),RES0(0.0),InnerIter(0),InnerItermax(500),
Literations(Liter), RESOR(4),EqnNames(4),SOR(4)
{
EqnNames[0] = "U-Eqn";
EqnNames[1] = "V-Eqn";
EqnNames[2] = "P-Eqn";
EqnNames[3] = "T-Eqn";

for(int i=0; i<SOR.size(); i++)
{
    SOR[i] = 0.2;
}

}

PCGSolver::~PCGSolver()
{
    //dtor
}

Fields::vectorField  PCGSolver::solve(Fields::vectorField& phi,Solution& sol_,FiniteMatrix::finiteMat& AS,FiniteMatrix::finiteMat& AN,FiniteMatrix::finiteMat& AW,
FiniteMatrix::finiteMat& AE,FiniteMatrix::finiteMat& AP,FiniteMatrix::finiteMat& SU, FiniteMatrix::finiteMat& SV,int& iterations,int& EqnNo)
{

            Fields::vectorField phitemp(phi.size(),vector<Fields>(phi[0].size()));

    int Liter = iterations;
//.....CALCULATE RESIDUAL AND OVERWRITE IT BY INTERMEDIATE VECTOR

        Residual=0.0;

        for(unsigned int i=1; i<phi.size()-1; i++) \
        {
          for(unsigned int j=1; j<phi[i].size()-1; j++)
            {
            RES[i][j].value  = SU[i][j].value - (AP[i][j].value*phi[i][j].value) -  (AE[i][j].value * phi[i+1][j].value) -
                        (AW[i][j].value * phi[i-1][j].value) - (AN[i][j].value*phi[i][j+1].value) -( AS[i][j].value*phi[i][j-1].value) ;
            Residual += abs(RES[i][j].value);
            }
        }
        //cout << endl;
        if(InnerIter==0)
        {
        RES0 = Residual;
        }

        //C.....PRECONDITIONING MATRIX DIAGONAL
        for(unsigned int i=1; i<phi.size()-1; i++) \
        {
         for(unsigned int j=1; j<phi[i].size()-1; j++)
          {
        D[i][j].value = 1.0/ (AP[i][j].value - (D[i-1][j].value * pow(AW[i][j].value,2.0))
                -  (D[i][j-1].value* pow(AS[i][j].value,2.0)));
          }
        }

        //C.....CALCULATION OF  ZK; INNER ITERATION LOOP

      double S0=1e20;

      for(int N=1; N<InnerItermax; N++)
      {

         //..FORWARD SUBSTITUTION
        for(unsigned int i=1; i<phi.size()-1; i++) \
        {
         for(unsigned int j=1; j<phi[i].size()-1; j++)
          {
         ZK[i][j].value = (RES[i][j].value - (AW[i][j].value * ZK[i-1][j].value)
            - (AS[i][j].value * ZK[i][j-1].value))*D[i][j].value;
          }
        }
        for(unsigned int i=1; i<phi.size()-1; i++) \
        {
         for(unsigned int j=1; j<phi[i].size()-1; j++)
          {
         ZK[i][j].value = ZK[i][j].value/(1.0e-20+D[i][j].value);
          }
        }
        // BACK SUBSTUTION
        double SK=0.0;
         for(unsigned int i=phi.size()-2; i >=1; --i)
        {
            for(unsigned int j=phi[i].size()-2; j >=1; --j)
            {
             ZK[i][j].value = (ZK[i][j].value - (AE[i][j].value * ZK[i+1][j].value)
                - (AN[i][j].value *ZK[i][j+1].value))*D[i][j].value;
                SK= SK+ RES[i][j].value *ZK[i][j].value;
            }
        }

        //C.....CALCULATE BETA AND NEW SEARCH VECTOR
       double Beta = SK/S0;
        for(unsigned int i=1; i<phi.size()-1; i++) \
        {
         for(unsigned int j=1; j<phi[i].size()-1; j++)
          {
         PK[i][j].value = ZK[i][j].value+Beta*PK[i][j].value;
          }
        }

    //C.....CALCULATE SCALAR PRODUCT (PK . A PK) AND ALPHA (A PK OVERWRITES ZK)
    double S2=0.0;
        for(unsigned int i=1; i<phi.size()-1; i++) \
        {
         for(unsigned int j=1; j<phi[i].size()-1; j++)
          {
         ZK[i][j].value = (AP[i][j].value*PK[i][j].value) + (AE[i][j].value*PK[i+1][j].value)
                + (AW[i][j].value*PK[i-1][j].value)  + (AN[i][j].value*PK[i][j+1].value) +  (AS[i][j].value*PK[i][j-1].value) ;
            S2=S2+ PK[i][j].value *ZK[i][j].value;
          }
        }
        double ALF= SK/S2;
  //C.....CALCULATE NEW RESIDUAL AND UPDATE VARIABLE
    Residual = 0.0;
        for(unsigned int i=1; i<phi.size()-1; i++) \
        {
         for(unsigned int j=1; j<phi[i].size()-1; j++)
          {
           phi[i][j].value = phi[i][j].value + ALF*PK[i][j].value;
           RES[i][j].value =   RES[i][j].value - ALF*ZK[i][j].value;
           Residual = Residual+ abs(RES[i][j].value);
          }
        }

        /*
            static bool runOnce= true;
            if(runOnce)
            {
            cout << EqnNames[EqnNo] << " with PCG It no: " << N << " Residual -> " << Residual << endl;
            }
        */

        S0=SK;

            //.CONVERGENCE CHECK (ON THE FINEST GRID ONLY)
        RSM=Residual/(RES0+1.0e-20);
        InnerIter++;

        double SORMAX=0.005;
        if(RSM < SORMAX)
        {
         cout << EqnNames[EqnNo] << " with PCG Solver Converged in " << InnerIter << " with Residual -> " <<  Residual << "and RSM " << RSM << endl;
         InnerIter=0;
         break;
        }


      } // end N=1:Liter


        return phitemp;
}


