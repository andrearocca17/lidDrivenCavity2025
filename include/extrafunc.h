#include <vector>
#include <iomanip>      // std::setprecision
#include <iostream>


void print2dvec(vector<vector<double> >& vecx)
{
    for(int i=0; i<vecx.size(); i++)
    {
    for(int j=0; j<vecx[0].size(); j++)
    {
                cout << std::setprecision(3) << vecx[i][j] << ' ';
    }
    cout << endl;
    }
}
