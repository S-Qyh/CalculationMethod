//
//  main.cpp
//  计算方法c++
//
//  Created by YH Q on 2022/3/25.
//

#include <iostream>
#include <math.h>
#include "iomanip"
#include "../EquationSolve.h"
using namespace::std;

double fun(double x){
    double y = x*x*x - sin(x) - 4*x + 1;
    return y;
}
double fun_1(double x){
    double y = x - fun(x) / ((3*x*x) - cos(x) - 4);
//    double y = pow(2*x + 5,1/3.0);

    return y;
}

int main(int argc, const char * argv[]) {
    std::cout << "Hello, World!\n";
    double result ,result_1;
    EquationSolve ls;
//二分法
//    cout<<"第一个解："<<endl;
//    result = ls.bisection(-5,0,0.001,fun);
//    cout<<endl;
//    cout<<"第二个解："<<endl;
//    result = ls.bisection(0,1,0.001,fun);
//    cout<<endl;
//    cout<<"第三个解："<<endl;
//    result = ls.bisection(1,2,0.001,fun);
//    cout<<endl;
    //迭代法
    cout<<"迭代法第一个解："<<endl;
    result_1 = ls.newtomIter(1,0.001,20,fun_1);
//
//////    cout<<fun(1.5);
//


//    result_1 = ls.newtomIter(1.5,0.0001,10,fun_1);

    return 0;
}
