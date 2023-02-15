//
// Created by YH Q on 2022/4/15.
//
#include "Forecast/Matrix.h"
#include "Interpolation.h"

using namespace std;
int main(){
    const double PI = 3.1415926535;
    Matrix<double> X(5,1),Y(5,1);
    X(1,1) = 0;
    X(2,1) = PI/6;
    X(3,1) = PI/4;
    X(4,1) = PI/3;
    X(5,1) = PI/2;

    Y(1,1) = 0;
    Y(2,1) = 0.5;
    Y(3,1) = 0.70711;
    Y(4,1) = 0.86603;
    Y(5,1) = 1;

    Interpolation le;
    double result;

//    le.LagrangeIter(X,Y,PI/5,result);


//    le.newtonIter_2(X,Y,PI/5,result);
    int flag = 1;
    le.newtonIter_2(X,Y,PI/5,result,flag);
    if(flag>0) {
        cout << "牛顿插值结果：" << endl;
        cout << result;
    } else{
        cout<<"牛顿插值法";
        cout<<"无解！"<<endl;
    }

    le.LagrangeIter(X,Y,PI/5,result);
    cout<<result<<endl;

    return 0;
}