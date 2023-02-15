//
// Created by YH Q on 2022/3/27.
//

#ifndef CLION_EQUATIONSOLVE_H
#define CLION_EQUATIONSOLVE_H

#include <iostream>
#include <math.h>
#include "iomanip"

using namespace std;

class EquationSolve {
private:
    double i = 0;
public:
    double bisection(double a, double b, double e, double (*pf)(double x));

    double newtomIter(double x0, double e, double N, double (*pf)(double x));

    void bisectionDispaly(double a, double b, int k) {

        cout << setprecision(6) << setiosflags(ios::left) << "K=" << k << "  [" << a << "," << b << "] e="
             << fabs(b - a) / 2 << "  mid = " << (b + a) / 2 <<endl;
    }

    void newtomIterDisplay(double x1, double x0, int k) {
        cout << setprecision(6) << setiosflags(ios::left) << "K=  " << k << "  Current Value = " << x1 << "  e = "
             << double(fabs(x1 - x0)) << endl;
    }
};


#endif //CLION_EQUATIONSOLVE_H
