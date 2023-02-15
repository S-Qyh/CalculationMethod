//
// Created by YH Q on 2022/4/15.
//
#include "Forecast/Matrix.h"
#ifndef CALCULATIONMETHOD_INTERPOLATION_H
#define CALCULATIONMETHOD_INTERPOLATION_H


class Interpolation {
public:
    Interpolation();

    void LagrangeBaseFun(Matrix<double> X,Matrix<double> Y,int k,double a,double &x1,int &flag);
    void LagrangeIter(Matrix<double> X,Matrix<double> Y,double a,double &y);
    void newtonIter_1(Matrix<double> X,Matrix<double> Y,double a,double &x);
    void newtonIter_2(Matrix<double> X,Matrix<double> Y,double a,double &x,int &flag);
};


#endif //CALCULATIONMETHOD_INTERPOLATION_H
