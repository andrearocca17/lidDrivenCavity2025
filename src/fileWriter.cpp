#include "fileWriter.h"
#include "Grid.h"
#include "Fields.h"

#include <cstdio>    // fopen, fprintf, fclose
#include <cstring>   // strncpy
#include <iostream>  // std::cerr
#include <string>    // std::string

fileWriter::fileWriter() {}
fileWriter::~fileWriter() {}

void fileWriter::writeUVP(std::string& name,
                          int time_,
                          Grid& Grid_,
                          Fields::vectorField& Utemp,
                          Fields::vectorField& Vtemp,
                          Fields::vectorField& Ptemp)
{
    // Build filename: name + time + ".dat"
    std::string timename = std::to_string(time_);
    std::string newfilename = name + timename + ".dat";

    // Convert to C‐string
    size_t bufSize = newfilename.size() + 1;
    char* cstr = new char[bufSize];
    std::strncpy(cstr, newfilename.c_str(), bufSize);
    cstr[bufSize - 1] = '\0';

    // Open file
    FILE* f = std::fopen(cstr, "w+t");
    if (!f) {
        std::cerr << "Failed to open file: " << newfilename << std::endl;
        delete[] cstr;
        return;
    }

    // Write header
    int NX = static_cast<int>(Utemp.size());
    int NY = static_cast<int>(Utemp[0].size());
    std::fprintf(f, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
    std::fprintf(f, "ZONE  F=POINT\n");
    std::fprintf(f, "I=%d, J=%d\n", NX, NY);

    // Write data
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            double x = Grid_.XC[i];
            double y = Grid_.YC[j];
            double u = Utemp[i][j].value;
            double v = Vtemp[i][j].value;
            double p = Ptemp[i][j].value;
            std::fprintf(f, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n",
                         x, y, u, v, p);
        }
    }

    // Clean up
    std::fclose(f);
    delete[] cstr;
}
