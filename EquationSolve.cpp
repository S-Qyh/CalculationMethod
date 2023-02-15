//
// Created by YH Q on 2022/3/27.
//

#include "EquationSolve.h"

#include <iostream>

using namespace std;


double EquationSolve::newtomIter(double x0, double e, double N, double (*pf)(double x)) {
    int k = 1;
    double x1;
    while (k < N) {

        x1 = (*pf)(x0);
        newtomIterDisplay(x1, x0, k);
        if (fabs(x1 - x0) < e) {
            cout << "方程的根为：" << x1;
            return x1;
        } else if (k < N) {
            k = k + 1;
            x0 = x1;
        }

    }
    return -999999.999;
}

double EquationSolve::bisection(double a, double b, double e, double (*pf)(double x)) {
    double x;
    if ((*pf)(a) * (*pf)(b) >= 0) {
        cout << "此区间无根" << endl;
        return NULL;
    }
    int i = 0;
    while (i < 1000) {
        x = (a + b) / 2;
        bisectionDispaly(a, b, i);
        if ((*pf)(x) == 0) {
            cout << "方程的解为：" << b << endl;
            return b;
        } else if ((*pf)(x) * (*pf)(a) >= 0) {
            a = x;
        } else if ((*pf)(x) * (*pf)(a) < 0) {
            b = x;
        }
        if (fabs(b - a) < e) {
            cout << "方程的解为：" << b;
            break;
        } else
            i = i + 1;
    }
    return b;
}