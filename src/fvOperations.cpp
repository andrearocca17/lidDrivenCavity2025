#include "fvOperations.h"

fvOperations::fvOperations(Grid& G, Fields& F):
Grid_(G),Field_(F)
{
 NI = Grid_.pNI();
 NJ = Grid_.pNJ();
 NIM = NI-1;
 NJM = NJ-1;

}
fvOperations::~fvOperations()
{
    //dtor
}
