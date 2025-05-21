#ifndef SIMPLE_H
#define SIMPLE_H
#include "Grid.h"
#include "Solution.h"
#include "FiniteMatrix.h"
#include "fileWriter.h"

class Equation;

class Simple
{
public:
    Simple(int& N_, int& M_, double& l_);
    virtual ~Simple();

    void assembleSolveMomentum();


private:
    void passGrid();
    void resolveU();
    void resolveV();
    void resolveP();

    int pressureGEast = 1;
    int pressureGNorth = 2;
    int solveIter = 1;
    double blendConvec = 1.0;
    Equation* UEqn = nullptr;
    Equation* VEqn = nullptr;
    Equation* PEqn = nullptr;
    fileWriter fileWriterr;
    string results = "OOPLid";
    //variabili per la definizione del problema specifico
  string EastB = "East"; string Northb = "North";
  string WestB = "West"; string Southb = "South";
  string timeMethod = "EULER2"; // either EULER2 or EULER3
  double ULID = 1.0;
  int NI;
  int NJ;
  Solution sol_;
  Fields FieldOper;
  Grid myGrid_;
  Fields::vectorField U;
  Fields::vectorField UO;
  Fields::vectorField UOO;
  Fields::vectorField V;
  Fields::vectorField VO;
  Fields::vectorField VOO;
  Fields::vectorField P;
  Fields::vectorField PP;
  Fields::vectorField DPX;
  Fields::vectorField DPY;
  Fields::vectorField F1;
  Fields::vectorField F2;
  FiniteMatrix::finiteMat AE;
  FiniteMatrix::finiteMat AW;
  FiniteMatrix::finiteMat AS;
  FiniteMatrix::finiteMat AN;
  FiniteMatrix::finiteMat AP;
  FiniteMatrix::finiteMat SU;
  FiniteMatrix::finiteMat SV;
  FiniteMatrix::finiteMat APU;
  FiniteMatrix::finiteMat APV;
  FiniteMatrix::finiteMat APUPtoE;
  FiniteMatrix::finiteMat APVPtoN;
  FiniteMatrix::finiteMat PAE;
  FiniteMatrix::finiteMat PAN;

};
#endif // SIMPLE_H