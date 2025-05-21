#include "FiniteVolume.h"
#include "forAllOperations.h"
#include "Solution.h"

namespace fvm {

FiniteMatrix::finiteMat dummyTerm(Fields::vectorField& vec)
{
    // create a matrix of the right size
    FiniteMatrix::finiteMat APtemp(
        vec.size(),
        std::vector<FiniteMatrix>(vec[0].size())
    );

    // towards East side
    forAllInternalUCVs(vec)
    {
        APtemp[i][j].aevalue   =  1.0 * vec[i][j].value;
        APtemp[i+1][j].awvalue = -1.0 * vec[i][j].value;
    }

    // towards North side
    forAllInternalVCVs(vec)
    {
        APtemp[i][j].anvalue   =  2.0 * vec[i][j].value;
        APtemp[i][j+1].asvalue = -2.0 * vec[i][j].value;
    }

    return APtemp;
}

FiniteMatrix::finiteMat diffusiveTerm(Fields::vectorField& vec)
{
    FiniteMatrix::finiteMat APtemp(
        vec.size(),
        std::vector<FiniteMatrix>(vec[0].size())
    );

    forAllInternalUCVs(vec)
    {
        APtemp[i][j].aevalue   = -(vec[i][j].visc * vec[i][j].Se) / vec[i][j].DXPtoE;
        APtemp[i+1][j].awvalue =  APtemp[i][j].aevalue;  // symmetry
    }

    forAllInternalVCVs(vec)
    {
        APtemp[i][j].anvalue   = -(vec[i][j].visc * vec[i][j].Sn) / vec[i][j].DYPtoN;
        APtemp[i][j+1].asvalue =  APtemp[i][j].anvalue;  // symmetry
    }

    return APtemp;
}

FiniteMatrix::finiteMat convDiffusive(
    Fields::vectorField& vec,
    Fields::vectorField& massFluxEast,
    Fields::vectorField& massFluxNorth,
    double& blendFac
) {
    FiniteMatrix::finiteMat APtemp(
        vec.size(),
        std::vector<FiniteMatrix>(vec[0].size())
    );

    // East face
    forAllInternalUCVs(vec)
    {
        double diff = (vec[i][j].visc * vec[i][j].Se) / vec[i][j].DXPtoE;
        double cE  = std::min(massFluxEast[i][j].value, 0.0);
        double cP  = std::max(massFluxEast[i][j].value, 0.0);

        double uds = cP * vec[i][j].value
                   + cE * vec[i+1][j].value;
        double cds = massFluxEast[i][j].value * (
                       vec[i+1][j].value * vec[i][j].FXE +
                       vec[i][j].value   * (1.0 - vec[i][j].FXE)
                     );
        double res = blendFac * (uds - cds);

        APtemp[i][j].aevalue   = cE - diff;
        APtemp[i+1][j].awvalue = -cP - diff;
        APtemp[i][j].svalue   +=  res;
        APtemp[i+1][j].svalue -=  res;
    }

    // North face
    forAllInternalVCVs(vec)
    {
        double diff = (vec[i][j].visc * vec[i][j].Sn) / vec[i][j].DYPtoN;
        double cN  = std::min(massFluxNorth[i][j].value, 0.0);
        double cP  = std::max(massFluxNorth[i][j].value, 0.0);

        double uds = cP * vec[i][j].value
                   + cN * vec[i][j+1].value;
        double cds = massFluxNorth[i][j].value * (
                       vec[i][j+1].value * vec[i][j].FYN +
                       vec[i][j].value   * (1.0 - vec[i][j].FYN)
                     );
        double res = blendFac * (uds - cds);

        APtemp[i][j].svalue   +=  res;
        APtemp[i][j+1].svalue -=  res;
    }

    return APtemp;
}

FiniteMatrix::finiteMat convectionTerm(
    Fields::vectorField& vec,
    Fields::vectorField& massFluxEast,
    Fields::vectorField& massFluxNorth,
    double& blendFac
) {
    // identical to convDiffusive but without diff terms
    FiniteMatrix::finiteMat APtemp(
        vec.size(),
        std::vector<FiniteMatrix>(vec[0].size())
    );

    forAllInternalUCVs(vec)
    {
        double cE  = std::min(massFluxEast[i][j].value, 0.0);
        double cP  = std::max(massFluxEast[i][j].value, 0.0);

        APtemp[i][j].aevalue   =  cE;
        APtemp[i+1][j].awvalue = -cP;
    }
    forAllInternalVCVs(vec)
    {
        double cN  = std::min(massFluxNorth[i][j].value, 0.0);
        double cP  = std::max(massFluxNorth[i][j].value, 0.0);

        APtemp[i][j].anvalue   =  cN;
        APtemp[i][j+1].asvalue = -cP;
    }

    return APtemp;
}

FiniteMatrix::finiteMat timeDerivative(
    Fields::vectorField& now,
    Fields::vectorField& old,
    Fields::vectorField& oldold,
    std::string& method
) {
    FiniteMatrix::finiteMat APtemp(
        now.size(),
        std::vector<FiniteMatrix>(now[0].size())
    );

    const double GAMT = (method == "EULER2" ? 1.0 : 0.5);  
    forAllInternalUCVs(now)
    {
        double valNow    =  now[i][j].value;
        double valOld    =  old[i][j].value;
        double valOldOld = oldold[i][j].value;

        double result = (method == "EULER3")
            ? (1.5*valNow - 2.0*valOld + 0.5*valOldOld)
            : (valNow - valOld);

        APtemp[i][j].value = result;
    }

    return APtemp;
}

FiniteMatrix::finiteMat pressureGrad(
    const Fields::vectorField& vec,
    int& direction_
) {
    FiniteMatrix::finiteMat APtemp(
        vec.size(),
        std::vector<FiniteMatrix>(vec[0].size())
    );

    forAllInternal(vec)
    {
        double DX = vec[i][j].X - vec[i-1][j].X;
        double DY = vec[i][j].Y - vec[i][j-1].Y;

        double e = (vec[i+1][j].value * vec[i][j].FXE)
                 + (vec[i][j].value * (1.0 - vec[i][j].FXE));
        double w = (vec[i][j].value   * vec[i-1][j].FXE)
                 + (vec[i-1][j].value * (1.0 - vec[i-1][j].FXE));
        double n = (vec[i][j+1].value * vec[i][j].FYN)
                 + (vec[i][j].value   * (1.0 - vec[i][j].FYN));
        double s = (vec[i][j].value   * vec[i][j-1].FYN)
                 + (vec[i][j-1].value * (1.0 - vec[i][j-1].FYN));

        double pE = (e - w)/DX;
        double pN = (n - s)/DY;

        if (direction_ == 1)
            APtemp[i][j].svalue = -pE * vec[i][j].volume;
        else
            APtemp[i][j].svalue = -pN * vec[i][j].volume;
    }

    return APtemp;
}

FiniteMatrix::finiteMat HTerm(
    const FiniteMatrix::finiteMat& APU1,
    const FiniteMatrix::finiteMat& APV1,
    Fields::vectorField& vec,
    Fields::vectorField& vec2
) {
    FiniteMatrix::finiteMat APtemp(
        APU1.size(),
        std::vector<FiniteMatrix>(APU1[0].size())
    );

    forAllInternalUCVs(APtemp)
    {
        APtemp[i][j].aevalue   = -vec[i][j].density * vec[i][j].Se * APU1[i][j].value * vec[i][j].Se;
        APtemp[i+1][j].awvalue =  APtemp[i][j].aevalue;
    }
    forAllInternalVCVs(APtemp)
    {
        APtemp[i][j].anvalue   = -vec2[i][j].density * vec2[i][j].Sn * APV1[i][j].value * vec[i][j].Sn;
        APtemp[i][j+1].asvalue =  APtemp[i][j].anvalue;
    }

    return APtemp;
}

FiniteMatrix::finiteMat divPhi(
    const Fields::vectorField& Feast,
    const Fields::vectorField& Fnorth
) {
    FiniteMatrix::finiteMat APtemp(
        Feast.size(),
        std::vector<FiniteMatrix>(Feast[0].size())
    );

    forAllInternal(APtemp)
    {
        APtemp[i][j].svalue =
            Feast[i-1][j].value - Feast[i][j].value
          + Fnorth[i][j-1].value - Fnorth[i][j].value;
    }

    return APtemp;
}

}  // namespace fvm
