#ifndef SOLUTION_H
#define SOLUTION_H


class Solution
{
    public:
        Solution();
        virtual ~Solution();

         double dt;
         double R;
         double visc;
         double density;
         int nsteps;
         int maxit;
         int Imonitor;
         int Jmonitor;
        double URFUvel;
        double URFVvel;
        double URFPressure;
         double URFU;
         double URFV;
         double URFP;
         double Alfa;



    protected:

    private:
};

#endif // SOLUTION_H
