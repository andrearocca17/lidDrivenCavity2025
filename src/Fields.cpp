#include "Fields.h"
#include <iomanip>           // std::setprecision
#include "forAllOperations.h"

Fields::Fields():
value(0.0)
{

}

Fields::Fields(int& NI_, int& NJ_) :
value(0.0), NI(NI_),NJ(NJ_), NIM(NI-1),NJM(NJ-1)
{
    //ctor
}

Fields::~Fields()
{
    //dtor
}

void Fields::getGridInfoPassed(Fields::vectorField& f,Grid& grid_, Solution& sol_)
{

    forAll(f)
    {
        f[i][j].X = grid_.X[i];
        f[i][j].Y = grid_.Y[j];
        f[i][j].XC = grid_.XC[i];
        f[i][j].YC = grid_.YC[j];
        f[i][j].FXE = grid_.FX[i];
        f[i][j].FYN = grid_.FY[j];
        f[i][j].FXP = 1.0-grid_.FX[i];
        f[i][j].FYP = 1.0-grid_.FY[j];
        f[i][j].visc =sol_.visc;
        f[i][j].density =sol_.density;
    }

    forAllInternal(f)
    {
            f[i][j].DXPtoE = grid_.XC[i+1] - grid_.XC[i];
            f[i][j].DYPtoN = grid_.YC[j+1] - grid_.YC[j];
            f[i][j].Se = (grid_.Y[j]-grid_.Y[j-1]);
            f[i][j].Sn = (grid_.X[i]-grid_.X[i-1]);
            f[i][j].volume = f[i][j].Se * f[i][j].Sn;
    }

}
void Fields::copyInternalField(Fields::vectorField& VecTo, Fields::vectorField& VecFrom)
{
        forAllInternal(VecTo)
        {
            VecTo[i][j].value =  VecFrom[i][j].value;

        }


}

void Fields::setVectorFieldGridFeatures(Fields::listVectorField& list1,Grid& g,Solution& s)
{
        for(int L=0; L<list1.size(); L++)
        {
            Fields::vectorField& tempvec(list1[L]);
            getGridInfoPassed(tempvec,g,s);
        }
}

void Fields::dirichletBoundary(Fields::vectorField& vec,string& direction,double& bvalue)
{
     if(direction=="East")
     {
            forEastBoundary(vec)
            {
                vec[i][j].value = bvalue;
            }
     }
     else if(direction=="West")
     {
            forWestBoundary(vec)
            {
                vec[i][j].value = bvalue;
            }
     }
     else if(direction=="South")
     {
            forSouthBoundary(vec)
            {
                vec[i][j].value = bvalue;
            }
     }
     else if(direction=="North")
    {
            forNorthBoundary(vec)
            {
                vec[i][j].value = bvalue;
            }
    }
    else
    {
        cout << " direction is wrong " << endl;
    }
}

void Fields::oscillatingValueBoundary(Fields::vectorField& vec,string& direction,double& omega, double& time,double& Velvalue)
{
     if(direction=="East")
     {
            forEastBoundary(vec)
            {
                vec[i][j].value = Velvalue*sin(omega*time);
            }
     }
     else if(direction=="West")
     {
            forWestBoundary(vec)
            {
              vec[i][j].value = Velvalue*sin(omega*time);
            }
     }
     else if(direction=="South")
     {
            forSouthBoundary(vec)
            {
               vec[i][j].value = Velvalue*sin(omega*time);
            }
     }
     else if(direction=="North")
    {
            forNorthBoundary(vec)
            {
                vec[i][j].value = Velvalue*sin(omega*time);
            }
    }
    else
    {
        cout << " direction is wrong " << endl;
    }
}


void Fields::extaraPolateZeroGrad(Fields::vectorField& vec,vector<double>& FXvec, vector<double>& FYvec,string& direction)
{

    // EAST
     if(direction=="East")
     {
          for(int j=1; j<vec[0].size()-1; j++)
            {
                int i=vec[0].size()-1;
                vec[i][j].value =   vec[i-1][j].value +  (vec[i-1][j].value- vec[i-2][j].value)*FXvec[i-1];
            }
     }

     else if(direction=="West")
     {
            for(int j=1; j<vec[0].size()-1; j++)
            {
                int i=0;
                vec[i][j].value =   vec[i+1][j].value +  (vec[i+1][j].value- vec[i+2][j].value)*FXvec[2];
            }
     }

     else if(direction=="South")
     {
         for(int i=1; i<vec.size()-1; i++)
            {
                int j=0;
                vec[i][j].value =   vec[i][j+1].value +  (vec[i][j+1].value- vec[i][j+2].value)*FYvec[2];
            }
     }

    else  if(direction=="North")
    {
         for(int i=1; i<vec.size()-1; i++)
            {
                int j=vec.size()-1;
                vec[i][j].value = vec[i][j-1].value +  (vec[i][j-1].value- vec[i][j-2].value)*(1.0-FYvec[j-1]);
            }
    }

    else
    {
        cout << " direction is wrong " << endl;
    }
}

void Fields::shiftSolutionInTime(Fields::vectorField& vecold, Fields::vectorField& vecnew)
{

    forAll(vecold)
    {
        vecold[i][j].value = vecnew[i][j].value;
    }


}

void Fields::initalizeFields(Fields::vectorField& vec, double& value)
{
    forAllInternal(vec)
    {
        vec[i][j].value = value;
    }
}


void Fields::computeCellCenterPressureGrad(Fields::vectorField& vec, Fields::vectorField& dvec,Fields::vectorField& dvec2)
{
forAllInternal(vec)
{
        double DX = vec[i][j].X-vec[i-1][j].X;
        double DY = vec[i][j].Y-vec[i][j-1].Y;


        double pressureEastFace = (vec[i+1][j].value*vec[i][j].FXE) + (vec[i][j].value* vec[i][j].FXP);
        double pressureWestFace = (vec[i][j].value*vec[i-1][j].FXE) + (vec[i-1][j].value*vec[i-1][j].FXP);
        double pressureNorthFace = (vec[i][j+1].value*vec[i][j].FYN) + (vec[i][j].value*vec[i][j].FYP);
        double pressureSouthFace = (vec[i][j].value*vec[i][j-1].FYN) + (vec[i][j-1].value*vec[i][j-1].FYP);
         dvec[i][j].value =  (pressureEastFace-pressureWestFace)/DX;
         dvec2[i][j].value = (pressureNorthFace-pressureSouthFace)/DY;

}
}


Fields::vectorField Fields::interpolatedFieldEast(Fields::vectorField& vec,Grid& myGrid_)
{
        Fields::vectorField temp(vec.size(),vector<Fields>(vec[0].size()));
        forAllInternalUCVs(vec)
        {
        double FXE = myGrid_.FX[i];
        double FXP = 1.0-FXE;
        temp[i][j].value = (vec[i+1][j].value*FXE) + (vec[i][j].value *FXP);
        }
        return temp;

}

Fields::vectorField Fields::interpolatedFieldNorth(Fields::vectorField& vec,Grid& myGrid_)
{
        Fields::vectorField temp(vec.size(),vector<Fields>(vec[0].size()));
        forAllInternalVCVs(vec)
        {
        double FYN = myGrid_.FY[j];
        double FYP = 1.0-FYN;
        temp[i][j].value = (vec[i][j+1].value*FYN) + (vec[i][j].value *FYP);
        }
        return temp;
}
Fields::vectorField Fields::cellFaceGradientEast(Fields::vectorField& vec,Grid& myGrid_)
{
        Fields::vectorField temp(vec.size(),vector<Fields>(vec[0].size()));
        forAllInternalUCVs(vec)
        {
        double DXPE = (myGrid_.XC[i+1]-myGrid_.XC[i]);
            temp[i][j].value =  ((vec[i+1][j].value) - (vec[i][j].value))/DXPE;
        }
        return temp;

}
Fields::vectorField Fields::cellFaceGradientNorth(Fields::vectorField& vec,Grid& myGrid_)
{
        Fields::vectorField temp(vec.size(),vector<Fields>(vec[0].size()));
      forAllInternalVCVs(vec)
        {
        double DYPN = (myGrid_.YC[j+1]-myGrid_.YC[j]);
            temp[i][j].value =  ((vec[i][j+1].value) - (vec[i][j].value))/DYPN;
        }
        return temp;

}

void Fields::computeEastMassFluxes(Fields::vectorField& vec,Fields::vectorField& corrU)
{
        forAllInternalUCVs(vec)
        {
            double sArea = vec[i][j].Se;
            double density = vec[i][j].density;
            vec[i][j].value = sArea*density*corrU[i][j].value;
        }

}
void Fields::computeNorthMassFluxes(Fields::vectorField& vec,Fields::vectorField& corrV)
{
      forAllInternalVCVs(vec)
        {
            double sArea = vec[i][j].Sn;
            double density = vec[i][j].density;

            vec[i][j].value = sArea*density*corrV[i][j].value;
        }
}


void Fields::print2dfield(Fields::vectorField& vec)
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
