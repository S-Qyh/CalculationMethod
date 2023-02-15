//
// Created by YH Q on 2022/4/11.
//

#ifndef CALCULATIONMETHOD_LINEAREQUATIONSSOLVE_H
#define CALCULATIONMETHOD_LINEAREQUATIONSSOLVE_H
#include "Matrix.h"

class LinearEquationsSolve {
public:
    //构造函数
    LinearEquationsSolve();
    //列高斯消元法
    void GaussCol(Matrix<double> &mat, Matrix<double> &b,int &flag,Matrix<double> &x);
    //选主元
    void ColPivot(Matrix<double> &A,Matrix<double> &B,int k);
    //高斯消元法
    void Gauss(Matrix<double> &A, Matrix<double> &B,Matrix<double> &x,int &flag);
    //LU分解
    void matrixLU(Matrix<double> &mat, Matrix<double> &b, int n);
    //雅可比迭代
    void jacobiIterator(Matrix<double> &A, Matrix<double> &B, double e, int M, int &flag, Matrix<double> X0, Matrix<double> &X);
    //高斯—塞的尔迭代法
    void gaussSD(Matrix<double> &A,Matrix<double> &B,double e,int M,int &flag,Matrix<double>X0,Matrix<double> &X);



};


#endif //CALCULATIONMETHOD_LINEAREQUATIONSSOLVE_H
