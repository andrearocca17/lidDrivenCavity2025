#include "SipSolSolver.h"
#include <iomanip>
#include <cmath>

SipSolSolver::SipSolSolver()
{

}
SipSolSolver::SipSolSolver(int& NI_, int& NJ_,int& Liter)
: NI(NI_),NJ(NJ_),NIM(NI-1),NJM(NJ-1),
UE(NI,Svector1d(NJ)),UN(NI,Svector1d(NJ)),LW(NI,Svector1d(NJ)),LS(NI,Svector1d(NJ)),LPR(NI,Svector1d(NJ)),
RES(NI,Svector1d(NJ)) , value(0.0),Residual(0.0),RSM(0.0),Literations(Liter), RESOR(4),EqnNames(4),SOR(4)
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

SipSolSolver::~SipSolSolver()
{
    //dtor

}

Fields::vectorField  SipSolSolver::solve(Fields::vectorField& phi,Solution& sol_,FiniteMatrix::finiteMat& AS,FiniteMatrix::finiteMat& AN,FiniteMatrix::finiteMat& AW,
FiniteMatrix::finiteMat& AE,FiniteMatrix::finiteMat& AP,FiniteMatrix::finiteMat& SU, FiniteMatrix::finiteMat& SV,int& iterations,int& EqnNo)
{

            Fields::vectorField phitemp(phi.size(),vector<Fields>(phi[0].size()));
//..COEFFICIENTS OF UPPER AND LOWER TRIANGULAR MATRICES
        for(unsigned int i=1; i<phi.size()-1; i++)
        {
                for(unsigned int j=1; j<phi[i].size()-1; j++)
          {
        LW[i][j].value = AW[i][j].value/(1.0+(sol_.Alfa*UN[i-1][j].value));
        LS[i][j].value = AS[i][j].value/(1.0+(sol_.Alfa*UE[i][j-1].value));
        double P1 = sol_.Alfa*LW[i][j].value * UN[i-1][j].value;
        double P2 = sol_.Alfa*LS[i][j].value * UE[i][j-1].value;
        LPR[i][j].value = 1.0/ (AP[i][j].value + P1 + P2 - LW[i][j].value*UE[i-1][j].value-LS[i][j].value*UN[i][j-1].value);
        UN[i][j].value = (AN[i][j].value - P1)*LPR[i][j].value;
        UE[i][j].value = (AE[i][j].value - P2)*LPR[i][j].value;

    }

    }


    int Liter = iterations;
//.....CALCULATE RESIDUAL AND OVERWRITE IT BY INTERMEDIATE VECTOR
    for(int L=0; L<Liter; L++)
    {
        Residual=0.0;

        for(unsigned int i=1; i<phi.size()-1; i++)
        {
                for(unsigned int j=1; j<phi[i].size()-1; j++)
          {
            RES[i][j].value  = SU[i][j].value - (AN[i][j].value*phi[i][j+1].value) -( AS[i][j].value*phi[i][j-1].value) -
                                (AE[i][j].value * phi[i+1][j].value) - (AW[i][j].value * phi[i-1][j].value)- (AP[i][j].value*phi[i][j].value);

            Residual += abs(RES[i][j].value);
            RES[i][j].value =   (RES[i][j].value - (LS[i][j].value* RES[i][j-1].value) - (LW[i][j].value* RES[i-1][j].value))*LPR[i][j].value;
                  //    cout << std::setprecision(3)  << RES[i][j].value << ' ';

        }
        //cout << endl;
        }

            double small=1e-20;
            if(L==0)
            {
            RESOR[EqnNo] = Residual;
            }
            RSM = Residual/(RESOR[EqnNo]+small);


      cout << EqnNames[EqnNo] << " Inner It: " << L << " and Residual --> " << Residual << " RSM " << RSM << endl;

    // Back Subsitution and Correction
    for(unsigned int i=phi.size()-2; i >=1; --i)
        {
     for(unsigned int j=phi[i].size()-2; j >=1; --j)
            {

            RES[i][j].value = RES[i][j].value - (UN[i][j].value*RES[i][j+1].value) - (UE[i][j].value*RES[i+1][j].value);
            phi[i][j].value = phi[i][j].value +  RES[i][j].value;
            }
        }

    }

        forAll(phitemp)
        {
            phitemp[i][j].value = phi[i][j].value;
        }
        return phitemp;
}


Fields::vectorField  SipSolSolver::solve(Fields::vectorField& phi,Solution& sol_,FiniteMatrix::finiteMat& AP,FiniteMatrix::finiteMat& SU,int& iterations,int& EqnNo)
{

            Fields::vectorField phitemp(phi.size(),vector<Fields>(phi[0].size()));
            FiniteMatrix::finiteMat AE(phi.size(),vector<FiniteMatrix>(phi[0].size()));
            FiniteMatrix::finiteMat AW(phi.size(),vector<FiniteMatrix>(phi[0].size()));
            FiniteMatrix::finiteMat AS(phi.size(),vector<FiniteMatrix>(phi[0].size()));
            FiniteMatrix::finiteMat AN(phi.size(),vector<FiniteMatrix>(phi[0].size()));

           for(unsigned int i=1; i<phi.size()-1; i++)
        {
                for(unsigned int j=1; j<phi[i].size()-1; j++)
          {
                AE[i][j].value = AP[i+1][j].value;
                AW[i][j].value = AP[i-1][j].value;
                AS[i][j].value = AP[i][j-1].value;
                AN[i][j].value = AP[i][j+1].value;

          }
      }



//..COEFFICIENTS OF UPPER AND LOWER TRIANGULAR MATRICES
        for(unsigned int i=1; i<phi.size()-1; i++)
        {
                for(unsigned int j=1; j<phi[i].size()-1; j++)
          {
        LW[i][j].value = AW[i][j].value/(1.0+(sol_.Alfa*UN[i-1][j].value));
        LS[i][j].value = AS[i][j].value/(1.0+(sol_.Alfa*UE[i][j-1].value));
        double P1 = sol_.Alfa*LW[i][j].value * UN[i-1][j].value;
        double P2 = sol_.Alfa*LS[i][j].value * UE[i][j-1].value;
        LPR[i][j].value = 1.0/ (AP[i][j].value + P1 + P2 - LW[i][j].value*UE[i-1][j].value-LS[i][j].value*UN[i][j-1].value);
        UN[i][j].value = (AN[i][j].value - P1)*LPR[i][j].value;
        UE[i][j].value = (AE[i][j].value - P2)*LPR[i][j].value;
    }
    }


    int Liter = iterations;
//.....CALCULATE RESIDUAL AND OVERWRITE IT BY INTERMEDIATE VECTOR
    for(int L=0; L<Liter; L++)
    {
        Residual=0.0;

        for(unsigned int i=1; i<phi.size()-1; i++)
        {
                for(unsigned int j=1; j<phi[i].size()-1; j++)
          {
            RES[i][j].value  = SU[i][j].value - (AN[i][j].value*phi[i][j+1].value) -( AS[i][j].value*phi[i][j-1].value) -
                                (AE[i][j].value * phi[i+1][j].value) - (AW[i][j].value * phi[i-1][j].value)- (AP[i][j].value*phi[i][j].value);

            Residual += abs(RES[i][j].value);
            RES[i][j].value =   (RES[i][j].value - (LS[i][j].value* RES[i][j-1].value) - (LW[i][j].value* RES[i-1][j].value))*LPR[i][j].value;
                  //    cout << std::setprecision(3)  << RES[i][j].value << ' ';

        }
        //cout << endl;
        }

            double small=1e-20;
            if(L==0)
            {
            RESOR[EqnNo] = Residual;
            }
            RSM = Residual/(RESOR[EqnNo]+small);


      cout << EqnNames[EqnNo] << " Inner It: " << L << " and Residual --> " << Residual << " RSM " << RSM << endl;

    // Back Subsitution and Correction
    for(unsigned int i=phi.size()-2; i >=1; --i)
        {
     for(unsigned int j=phi[i].size()-2; j >=1; --j)
            {

            RES[i][j].value = RES[i][j].value - (UN[i][j].value*RES[i][j+1].value) - (UE[i][j].value*RES[i+1][j].value);
            phi[i][j].value = phi[i][j].value +  RES[i][j].value;
            }
        }

    }

        forAll(phitemp)
        {
            phitemp[i][j].value = phi[i][j].value;
        }
        return phitemp;


}

void SipSolSolver::printSvector(SipSolSolver::Svector& vec)
{
    for(int i=0; i<vec.size(); i++)
    {
    for(int j=0; j<vec[i].size(); j++)
    {
    cout << std::setprecision(3) << vec[i][j].value << ' ';
    }
    cout << endl;
    }

}
