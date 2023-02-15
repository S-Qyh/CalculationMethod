//
// Created by YH Q on 2022/5/4.
//

#include <cmath>
#include "NumericalIntegration.h"
#include "iostream"

using namespace std;

NumericalIntegration::NumericalIntegration(double (*p)(double x, double y), double (*Base)(double x)) {
    this->f = p;
    this->Base = Base;
}
NumericalIntegration::NumericalIntegration(double (*p)(double)) {
    this->Base = p;
}

void NumericalIntegration::trapezoidal(double a, double b, double (*pf)(double), double &y) {
    y = (1.0 / 2) * (b - a) * ((*pf)(a) + (*pf)(b));
}

void NumericalIntegration::rectangular(double a, double b, double (*pf)(double), double &y) {
    y = (b - a) * (*pf)((a + b) / 2);
}

void NumericalIntegration::simpson(double a, double b, double (*pf)(double), double &y) {
    y = 1.0 / 6 * (b - a) * ((*pf)(a) + 4 * (*pf)((a + b) / 2) + (*pf)(b));
}

void NumericalIntegration::compoundTrapezoid(double a, double b, int N, double (*pf)(double), double &y) {
    double h = (b - a) / N;
    y = 0;
    for (int k = 1; k < N; k++) {
        y = y + (*pf)(a + k * h);
    }
    y = h * ((*pf)(a) + 2 * y + (*pf)(b)) / 2;
}

void NumericalIntegration::compoundSimpson(double a, double b, int N, double (*pf)(double), double &y) {
    double h = (b - a) / N;
    double x = a + h / 2;
    double S1 = (*pf)(x);
    double S2 = 0;
    for (int k = 1; k < N; k++) {
        S1 = S1 + (*pf)(a + k * h + h / 2);
        S2 = S2 + (*pf)(a + k * h);
    }
    y = h * ((*pf)(a) + 4 * S1 + 2 * S2 + (*pf)(b)) / 6;
}

void NumericalIntegration::varStepTrapezoid(double a, double b, double e, double (*pf)(double), double &y) {
    double h = b - a;
    double T1 = h / 2 * ((*pf)(a) + (*pf)(b));
    while (true) {
        double S = 0;
        double x = a + h / 2;
        do {
            S = S + (*pf)(x);
            x = x + h;
        } while (x < b);
        double T2 = (T1 + h * S) / 2;
        if (fabs(T2 - T1) <= e) {
            y = T2;
            break;
        } else {
            T1 = T2;
            h = h / 2;
        }
    }
}

void NumericalIntegration::fourRungeKutta(double x0, double y0, double h, double N, double (*pf)(double, double),
                                          double &y, double (*funBase)(double x)) {
    double n = 1;
    NumericalIntegration::display(x0, y0, funBase);
    while (true) {
        double x1 = x0 + h;
        double k1 = (*pf)(x0, y0);
        double k2 = (*pf)(x0 + h / 2, y0 + h * k1 / 2);
        double k3 = (*pf)(x0 + h / 2, y0 + h * k2 / 2);
        double k4 = (*pf)(x1, y0 + h * k3);
        double y1 = y0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        NumericalIntegration::display(x1, y1, funBase);
        if (n == N) {
            y = y1;
            break;
        } else {
            n = n + 1;
            x0 = x1;
            y0 = y1;
        }
    }
}

void
NumericalIntegration::imporveEuler(double x0, double y0, double h, double N, double (*pf)(double, double), double &y,
                                   double (*funBase)(double x)) {
    double n = 1;
    NumericalIntegration::display(x0, y0, funBase);
    while (true) {
        double x1 = x0 + h;
        double yp = y0 + h * (*pf)(x0, y0);
        double yc = y0 + h * (*pf)(x1, yp);
        double y1 = (yp + yc) / 2;
        NumericalIntegration::display(x1, y1, funBase);
        if (n == N) {
            y = y1;
            break;
        } else {
            n = n + 1;
            x0 = x1;
            y0 = y1;
        }
    }

}

void
NumericalIntegration::imporveEuler(double x0, double y0, double h, double N, double &y) {
    double n = 1;
    NumericalIntegration::display(x0, y0, this->Base);
    while (true) {
        double x1 = x0 + h;
        double yp = y0 + h * (*f)(x0, y0);
        double yc = y0 + h * (*f)(x1, yp);
        double y1 = (yp + yc) / 2;
        NumericalIntegration::display(x1, y1, Base);
        if (n == N) {
            y = y1;
            break;
        } else {
            n = n + 1;
            x0 = x1;
            y0 = y1;
        }
    }

}


void NumericalIntegration::fourRungeKutta(double x0, double y0, double h, double N, double &y) {
    double n = 1;
    NumericalIntegration::display(x0, y0);
    while (true) {
        double x1 = x0 + h;
        double yp = y0 + h * (*f)(x0, y0);
        double yc = y0 + h * (*f)(x1, yp);
        double y1 = (yp + yc) / 2;
        NumericalIntegration::display(x1, y1);
        if (n == N) {
            y = y1;
            break;
        } else {
            n = n + 1;
            x0 = x1;
            y0 = y1;
        }
    }
}

void NumericalIntegration::display(double x, double y, double (*pf)(double)) {
    cout << "x = " << x << "\t\t" << "y = " << y << "\t\t" << "y' = " << pf(x) << "\t\t" << "|y'-y| = "
         << fabs(pf(x) - y) << endl;
}

void NumericalIntegration::display(double x, double y) {
    cout << "x = " << x << "\t\t" << "y = " << y << "\t\t" << "y' = " << this->Base(x) << "\t\t" << "|y'-y| = "
         << fabs(this->Base(x) - y) << endl;
}

////构造的时候传入函数，这里不需要传入
void NumericalIntegration::compoundTrapezoid(double a, double b, int N, double &y) {
    double h = (b - a) / N;
    y = 0;
    for (int k = 1; k < N; k++) {
        y = y + this->Base(a + k * h);
    }
    y = h * (this->Base(a) + 2 * y + this->Base(b)) / 2;
}
////构造的时候传入函数，这里不需要传入
void NumericalIntegration::varStepTrapezoid(double a, double b, double e, double &y) {
    double h = b - a;
    int n = 1;
    double T1 = h / 2 * (this->Base(a) + this->Base(b));
    cout<<"n = "<<n<<"\tI = "<<T1<<endl;
    while (true) {
        double S = 0;
        double x = a + h / 2;
        do {
            S = S + this->Base(x);
            x = x + h;
        } while (x < b);
        n++;
        double T2 = (T1 + h * S) / 2;
        cout<<"n = "<<n<<"\tI = "<<T2<<endl;
        if (fabs(T2 - T1) <= e) {
            y = T2;
            break;
        } else {
            T1 = T2;
            h = h / 2;
        }
    }
}
////构造的时候传入函数，这里不需要传入
void NumericalIntegration::compoundSimpson(double a, double b, int N, double &y) {
    double h = (b - a) / N;
    double x = a + h / 2;
    double S1 = this->Base(x);
    double S2 = 0;
    for (int k = 1; k < N; k++) {
        S1 = S1 + this->Base(a + k * h + h / 2);
        S2 = S2 + this->Base(a + k * h);
    }
    y = h * (this->Base(a) + 4 * S1 + 2 * S2 + this->Base(b)) / 6;
}
////重载梯形积分公式，不传函数，构造的时候传入
void NumericalIntegration::trapezoidal(double a, double b, double &y) {
    y = (1.0 / 2) * (b - a) * (this->Base(a) + this->Base(b));
}
////重载中矩形公式，不传函数，构造传入
void NumericalIntegration::rectangular(double a, double b, double &y) {
    y = (b - a) * this->Base((a + b) / 2);
}
////重载辛普森，不传入函数。构造传入
void NumericalIntegration::simpson(double a, double b, double &y) {
    y = 1.0 / 6 * (b - a) * (this->Base(a) + 4 * this->Base((a + b) / 2) + this->Base(b));
}








