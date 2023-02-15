//
// Created by YH Q on 2022/5/4.
//
#include "NumericalIntegration.h"
#include "iostream"
#include <cmath>

using namespace std;

double fun(double x, double y) {
    return -10*y + x;
}

double funBase(double x) {
    return cos(x) + sin(x);
}

double fun1(double x) {
    return 1.0 / (1 + x * x);
}

double fun2(double x) {
    return exp(-x);
}

double fun3(double x){
    return sin(x)/x;
}

int main() {
    NumericalIntegration q;
    double y;
    cout<<"步长为：0.1 时候："<<endl;
    q.imporveEuler(0, 0, 0.1, 4, fun, y, funBase);
    cout << y << endl;
    cout << endl;
    cout<<"步长为：0.4 时候："<<endl;
    q.imporveEuler(0, 1, 0.4, 5, fun, y, funBase);
    cout << y << endl;




    cout<<"四阶龙哥库塔步长为0。2时候：" <<endl;
    q.fourRungeKutta(0, 0, 0.2, 2, fun, y, funBase);
    cout << y << endl;
    q.fourRungeKutta(0, 1, 0.4, 5, fun, y, funBase);
    cout << y << endl;

    NumericalIntegration s(fun, funBase);
    s.imporveEuler(0, 1, 0.2, 10, y);
    cout << y << endl;

    s.fourRungeKutta(0, 1, 0.4, 5, y);
    cout << y << endl;
    cout << endl;

    NumericalIntegration m(fun1);
    cout<<"变步长梯形求解："<<endl;
    m.varStepTrapezoid(0, 1, 0.0001, y);
    cout << y << endl;

    m.compoundSimpson(0, 1, 30, y);
    cout << y << endl;

    m.compoundSimpson(0, 1, 50, y);
    cout << y << endl;

    m.compoundSimpson(0, 1, 70, y);
    cout << y << endl;

    m.compoundSimpson(0, 1, 100, y);
    cout << y << endl;

    cout << "复化梯形：" << endl;
    m.compoundTrapezoid(0, 1, 30, y);
    cout << y << endl;

    m.compoundTrapezoid(0, 1, 50, y);
    cout << y << endl;

    m.compoundTrapezoid(0, 1, 70, y);
    cout << y << endl;

    m.compoundTrapezoid(0, 1, 100, y);
    cout << y << endl;

    m.compoundTrapezoid(0, 1, 150, y);
    cout << y << endl;

    NumericalIntegration t(fun2);
    t.trapezoidal(0, 1, fun2, y);
    cout << y << endl;
    t.simpson(0, 1, fun2, y);
    cout << y << endl;

    t.compoundTrapezoid(0, 1, 10, y);
    cout << y << endl;

    t.compoundSimpson(0, 1, 10, y);
    cout << y << endl;

    t.varStepTrapezoid(0, 1, 0.001, y);
    cout << y << endl;

    NumericalIntegration l;
    cout<<"start"<<endl;
    l.varStepTrapezoid(0,1,0.0001,y);
    cout<<y<<endl;


    return 0;
}