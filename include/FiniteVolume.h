#include "Grid.h"
#include "Fields.h"
#include "FiniteMatrix.h"

namespace fvm
{
    FiniteMatrix::finiteMat dummyTerm(Fields::vectorField& vec);
    FiniteMatrix::finiteMat diffusiveTerm(Fields::vectorField& vec);
    FiniteMatrix::finiteMat convDiffusive(Fields::vectorField& vec, Fields::vectorField& massFluxEast, Fields::vectorField& massFluxNorth, double& blendFac);
    FiniteMatrix::finiteMat convectionTerm(Fields::vectorField& vec, Fields::vectorField& massFluxEast, Fields::vectorField& massFluxNorth, double& blendFac);
    FiniteMatrix::finiteMat timeDerivative(Fields::vectorField& VecNow, Fields::vectorField& VecOld, Fields::vectorField& VecOldOld, string& method);
    FiniteMatrix::finiteMat pressureGrad(const Fields::vectorField& vec, int& direction_);
    FiniteMatrix::finiteMat HTerm(const FiniteMatrix::finiteMat& APU1, const FiniteMatrix::finiteMat& APV1, Fields::vectorField& vec, Fields::vectorField& vec2);
    FiniteMatrix::finiteMat divPhi(const Fields::vectorField& Feast, const Fields::vectorField& Fnorth);
}   //end namespace
