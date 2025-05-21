#include "Solution.h"

Solution::Solution():
dt(0.001),R(1.0),visc(1e-03),density(1.0),
nsteps(1),maxit(2),URFUvel(0.8),URFVvel(0.8),URFPressure(0.2),Imonitor(1),Jmonitor(1)
,Alfa(0.92)
{
    URFU = 1.0/URFUvel;
    URFV = 1.0/URFVvel;
    URFP = 1.0/URFPressure;

    //ctor
}

Solution::~Solution()
{
    //dtor
}
