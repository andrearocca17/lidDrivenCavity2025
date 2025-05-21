#include "FiniteMatrix.h"
#include "forAllOperations.h"
#include <iomanip>      // std::setprecision
FiniteMatrix::FiniteMatrix() : value(0.0) ,aevalue(0.0),awvalue(0.0),asvalue(0.0),anvalue(0.0), svalue(0.0)
{
    //ctor
}

FiniteMatrix::~FiniteMatrix()
{
    //dtor
}

FiniteMatrix::finiteMat FiniteMatrix::interpolatedFieldEast(FiniteMatrix::finiteMat& vec,Grid& myGrid_)
{
        FiniteMatrix::finiteMat temp(vec.size(),vector<FiniteMatrix>(vec[0].size()));
        forAllInternalUCVs(vec)
        {
        double FXE = myGrid_.FX[i];
        double FXP = 1.0-FXE;
        temp[i][j].value =  (vec[i+1][j].value*FXE) + (vec[i][j].value *FXP);
        }
        return temp;

}

FiniteMatrix::finiteMat FiniteMatrix::interpolatedFieldNorth(FiniteMatrix::finiteMat& vec,Grid& myGrid_)
{
        FiniteMatrix::finiteMat temp(vec.size(),vector<FiniteMatrix>(vec[0].size()));

        forAllInternalVCVs(vec)
        {
        double FYN = myGrid_.FY[j];
        double FYP = 1.0-FYN;
        temp[i][j].value =  (vec[i][j+1].value*FYN) + (vec[i][j].value *FYP);
        }
        return temp;
}
Fields::vectorField FiniteMatrix::correctedFaceVelocityEast(Fields::vectorField& interpolatedCellFaceCenterVelocity,Fields::vectorField& cellFacePressureGrad, Fields::vectorField& DPXfield,
        FiniteMatrix::finiteMat& APinterpolated,Grid& myGrid_)
{
       Fields::vectorField temp(interpolatedCellFaceCenterVelocity.size(),vector<Fields>(interpolatedCellFaceCenterVelocity[0].size()));

        forAllInternalUCVs(temp)
        {
            double DXPE = (myGrid_.XC[i+1]-myGrid_.XC[i]);
            double sArea =  (myGrid_.Y[j]-myGrid_.Y[j-1]);
            double volume =   DXPE*sArea;
        temp[i][j].value = interpolatedCellFaceCenterVelocity[i][j].value - (APinterpolated[i][j].value*volume* (cellFacePressureGrad[i][j].value -DPXfield[i][j].value) );

        }
    return temp;
}

Fields::vectorField FiniteMatrix::correctedFaceVelocityNorth(Fields::vectorField& interpolatedCellFaceCenterVelocity,Fields::vectorField& cellFacePressureGrad, Fields::vectorField& DPYfield,
        FiniteMatrix::finiteMat& APinterpolated, Grid& myGrid_)
{
       Fields::vectorField temp(interpolatedCellFaceCenterVelocity.size(),vector<Fields>(interpolatedCellFaceCenterVelocity[0].size()));

        forAllInternalVCVs(temp)
        {

        double DYPN = (myGrid_.YC[j+1]-myGrid_.YC[j]);
        double sArea = (myGrid_.X[i]-myGrid_.X[i-1]);
        double volume =  DYPN*sArea;
        temp[i][j].value = interpolatedCellFaceCenterVelocity[i][j].value - (APinterpolated[i][j].value*volume * (cellFacePressureGrad[i][j].value -DPYfield[i][j].value) );

        }
    return temp;
}
void FiniteMatrix::correctEastMassFluxes(Fields::vectorField& F1,Fields::vectorField& PP, FiniteMatrix::finiteMat& AE)
{
    //10.1 East
    for(int i=1; i< F1.size()-2; i++)
     {
        for(int j=1; j< F1[0].size()-1; j++)
        {
        F1[i][j].value = F1[i][j].value + (AE[i][j].value*(PP[i+1][j].value -PP[i][j].value));
        }
     }
}

void FiniteMatrix::correctNorthMassFluxes(Fields::vectorField& F2,Fields::vectorField& PP,FiniteMatrix::finiteMat& AN)
{
    //10.2 North
    for(int i=1; i< F2.size()-1; i++)
     {
        for(int j=1; j< F2[0].size()-2; j++)
      {
        F2[i][j].value = F2[i][j].value + (AN[i][j].value*(PP[i][j+1].value -PP[i][j].value));
       }
     }
}

void FiniteMatrix::print2dfield(FiniteMatrix::finiteMat& vec)
{
    for(int i=0; i<vec.size(); i++)
    {
    for(int j=0; j<vec[i].size(); j++)
    {
     std::cout << std::setprecision(3) << vec[i][j].value << ' ';
    }
    std::cout << endl;
    }

}
void FiniteMatrix::print2dfieldsource(FiniteMatrix::finiteMat& vec)
{
    for(int i=0; i<vec.size(); i++)
    {
    for(int j=0; j<vec[i].size(); j++)
    {
     std::cout << std::setprecision(3) << vec[i][j].svalue << ' ';
    }
    std::cout << endl;
    }

}




// OVERLOADED OPERATORS




FiniteMatrix::finiteMat operator+=(const FiniteMatrix::finiteMat& A, FiniteMatrix::finiteMat& B)
{
 forAllInternal(B)
 {
  B[i][j].value =   A[i][j].value +   B[i][j].value;

 }
 return B;
}

FiniteMatrix::finiteMat operator-=(const FiniteMatrix::finiteMat& A,  FiniteMatrix::finiteMat& B)
{
 forAllInternal(B)
 {
  B[i][j].value =   A[i][j].value - B[i][j].value;
 }
 return  B;
}

FiniteMatrix::finiteMat operator+(const FiniteMatrix::finiteMat& A, const FiniteMatrix::finiteMat& B)
{
 FiniteMatrix::finiteMat result(A);
 forAllInternal(A)
 {
  result[i][j].value += B[i][j].value;
  result[i][j].aevalue += B[i][j].aevalue;
  result[i][j].anvalue += B[i][j].anvalue;
  result[i][j].asvalue += B[i][j].asvalue;
  result[i][j].awvalue += B[i][j].awvalue;
  result[i][j].svalue += B[i][j].svalue;
 }
 return  result;
}

FiniteMatrix::finiteMat operator-(const FiniteMatrix::finiteMat& A, const FiniteMatrix::finiteMat& B)
{
 FiniteMatrix::finiteMat result(A);
 forAllInternal(A)
 {
  result[i][j].value -= B[i][j].value;
  result[i][j].aevalue -= B[i][j].aevalue;
  result[i][j].anvalue -= B[i][j].anvalue;
  result[i][j].asvalue -= B[i][j].asvalue;
  result[i][j].awvalue -= B[i][j].awvalue;
  result[i][j].svalue -= B[i][j].svalue;
 }
 return  result;
}

FiniteMatrix::finiteMat operator&&(const FiniteMatrix::finiteMat& A, const FiniteMatrix::finiteMat& B)
{
 FiniteMatrix::finiteMat result(A);
 forAllInternal(A)
 {
  result[i][j].value *= B[i][j].value;
 }
 return  result;
}
FiniteMatrix::finiteMat operator&&(const FiniteMatrix::finiteMat& A, const Fields::vectorField& B)
{
 FiniteMatrix::finiteMat result(A);
 forAllInternal(A)
 {
  result[i][j].value *= B[i][j].value;
 }
 return  result;
}
FiniteMatrix::finiteMat operator*(const double& dble, const FiniteMatrix::finiteMat& B)
{
 FiniteMatrix::finiteMat result(B);
 forAllInternal(B)
 {
  result[i][j].value =dble* B[i][j].value;
 }
 return  result;
}

FiniteMatrix::finiteMat operator*(const FiniteMatrix::finiteMat& B,const double& dble)
{
 FiniteMatrix::finiteMat result(B);
 forAllInternal(B)
 {
  result[i][j].value =dble* B[i][j].value;
 }
 return  result;
}


FiniteMatrix::finiteMat operator/(const double& dble, const FiniteMatrix::finiteMat& B)
{
 FiniteMatrix::finiteMat result(B);
 forAllInternal(B)
 {
  result[i][j].value =dble/B[i][j].value;
 }
 return  result;
}
FiniteMatrix::finiteMat operator/(const FiniteMatrix::finiteMat& B,const double& dble)
{
 FiniteMatrix::finiteMat result(B);
 forAllInternal(B)
 {
  result[i][j].value =B[i][j].value/dble;
 }
 return  result;
}
