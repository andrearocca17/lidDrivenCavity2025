#ifndef FILEWRITER_H
#define FILEWRITER_H
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "Grid.h"
#include "Fields.h"
#include <cstring>
#include <sstream>
using namespace std;
using std::vector;
using std::string;



class fileWriter
{
public:
    fileWriter();

    virtual ~fileWriter();

      void writeUVP(string&,int,Grid&, Fields::vectorField&, Fields::vectorField&,Fields::vectorField&);


};

#endif // FILEWRITER_H
