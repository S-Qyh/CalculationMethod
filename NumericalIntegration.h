//
// Created by YH Q on 2022/5/4.
//

#ifndef CALCULATIONMETHOD_NUMERICALINTEGRATION_H
#define CALCULATIONMETHOD_NUMERICALINTEGRATION_H


class NumericalIntegration {
public:
    NumericalIntegration(double (*p)(double x, double y), double (*Base)(double x));
    NumericalIntegration(double (*p)(double x));

    NumericalIntegration() {};

    void trapezoidal(double a, double b, double (*pf)(double x), double &y);  //梯形公式
    void trapezoidal(double a,double b,double &y);   // 梯形公式不传参数

    void rectangular(double a, double b, double (*pf)(double x), double &y);  //中矩形公式
    void rectangular(double a, double b, double &y);  //中矩形公式不传函数

    void simpson(double a, double b, double (*pf)(double x), double &y);      //辛普森公式
    void simpson(double a, double b, double &y);      //辛普森公式不传函数

    void compoundTrapezoid(double a, double b, int N, double (*pf)(double x), double &y);  //复化梯形
    void compoundTrapezoid(double a,double b,int N,double &y); //复化梯形不传函数

    void compoundSimpson(double a, double b, int N, double (*pf)(double x), double &y);    //复化辛普森
    void compoundSimpson(double a, double b, int N,double &y);    //复化辛普森不传函数

    void varStepTrapezoid(double a, double b, double e, double (*pf)(double x), double &y);    //变步长梯形公式
    void varStepTrapezoid(double a, double b, double e,double &y);   //变步长梯形公式不传函数

    void
    fourRungeKutta(double x0, double y0, double h, double N, double (*pf)(double x, double y), double &y,
                   double (*funBase)(double x));  // 四阶龙哥库塔
    void
    fourRungeKutta(double x0, double y0, double h, double N, double &y);  // 四阶龙哥库塔不加函数


    void imporveEuler(double x0, double y0, double h, double N, double (*pf)(double x, double y), double &y,
                      double (*funBase)(double x));     //欧拉法
    void imporveEuler(double x0, double y0, double h, double N, double &y);


    void display(double x, double y, double (*pf)(double x)); //显示

    void display(double x, double y); //显示不加函数


private:
    double (*f)(double x, double y);

    double (*Base)(double x);


};


#endif //CALCULATIONMETHOD_NUMERICALINTEGRATION_H
