#ifndef MATRIXCLASS_H
#define MATRIXCLASS_H

#include <vector>
#include <iostream>
#include "forAllOperations.h"
using namespace std;

#include "Fields.h"
using std::vector;


class MatrixClass
{
    public:
        MatrixClass();
        virtual ~MatrixClass();

            typedef vector<vector<MatrixClass> > mat2d;

            double value;

            void operator+(const MatrixClass&);
            void operator-(const MatrixClass&);
       //     void operator==(const MatrixClass&);

            // global operators
       friend MatrixClass::mat2d operator==(MatrixClass::mat2d&, const MatrixClass::mat2d&);
    protected:

    private:
};

#endif // MATRIXCLASS_H
