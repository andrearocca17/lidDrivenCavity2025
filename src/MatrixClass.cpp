#include "MatrixClass.h"


MatrixClass::MatrixClass():value(0.0)
{
    //ctor
}

MatrixClass::~MatrixClass()
{
    //dtor
}

void MatrixClass::operator+(const MatrixClass& A)
{
       value += A.value;
}

void MatrixClass::operator-(const MatrixClass& A)
{
       value -= A.value;
}


MatrixClass::mat2d operator==(MatrixClass::mat2d& A, const MatrixClass::mat2d& B)
{
 forAllInternal(A)
 {
  A[i][j].value =   A[i][j].value +   B[i][j].value;
 }
 return  A;
}
