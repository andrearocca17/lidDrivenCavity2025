#ifndef FVOPERATIONS_H
#define FVOPERATIONS_H

#include "Fields.h"
#include "Grid.h"


class fvOperations
{
    public:
        fvOperations(Grid&, Fields&);
        virtual ~fvOperations();

        int NI,NJ,NIM,NJM;


    protected:

    private:
        Grid Grid_;
        Fields Field_;

};

#endif // FVOPERATIONS_H
