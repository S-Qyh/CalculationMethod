//
// Created by YH Q on 2022/4/15.
//

#include "Interpolation.h"
#include "set"

Interpolation::Interpolation() {}


void Interpolation::LagrangeBaseFun(Matrix<double> X, Matrix<double> Y, int k, double a,double &x1,int &flag) {
    set<double> myset;
    for(int m = 1;m<=X.getNumRow();m++){
        myset.insert(X(m,1));
    }
    if(myset.size()!=X.getNumRow()){
        cout<<"不符合求值条件"<<endl;
        flag = 0;
        return;
    }
    double prod = 1;
    int n = X.getNumRow();
    for(int i = 1;i<=n;i++){
        if (i != k){
            prod = prod*(((a-X(i,1)))/ (X(k,1) - X(i,1)));
        }
    }
    cout<<"当前点"<<"  "<<k<<"  "<<"拉格朗日差值基函数值为：\n";
    cout<<prod<<endl;
    x1 = prod;
}

void Interpolation::LagrangeIter(Matrix<double> X, Matrix<double> Y, double a,double &y) {
    int n = X.getNumRow();
    double sum = 0;
    int flag;
    for (int i = 1; i <= n; i++) {
        double x1;
        LagrangeBaseFun(X,Y,i,a,x1,flag);
        if(flag>0) {
            sum += x1 * Y(i, 1);
        } else{
            cout<<"拉格朗日求法";
            cout<<"无解！"<<endl;
            return;
        }
    }
    y = sum;
    cout<<"拉格朗日结果："<<endl;
    cout << y << endl;
    cout<<endl;
}

////普通方法
void Interpolation::newtonIter_1(Matrix<double> X, Matrix<double> Y, double a,double &x) {
    double r = 0;
    int row = X.getNumRow();
    Matrix<double> F(row,row);
    F = Y;
    for(int j = 2;j<=row;j++){
        for(int i=j;i<=row;i++){
            F(i,j) = (F(i,j-1) - F(i-1,j-1))/(X(i,1) - X(i-j+1,1));
        }
    }
    Matrix<double>P(row,1);
    for(int i = 1;i<= row;i++){
        P(i,1) = F(i,i);
    }
    cout<<"差商表："<<endl;
    P.printMatrix(P);
    for(int i = 1;i<=row;i++){
        double temp = 1;
        for(int j = 1;j<i;j++){
            temp *= (a-X(j,1));
        }
        temp = temp * P(i,1);
        r += temp;
    }
    x = r;
}
////节省内存空间写法，从后往前计算差商表。
void Interpolation::newtonIter_2(Matrix<double> X, Matrix<double> Y, double a,double &x,int &flag) {
    set<double> myset;
    for(int m = 1;m<=X.getNumRow();m++){
        myset.insert(X(m,1));
    }
    if(myset.size()!=X.getNumRow()){
        cout<<"不符合求值条件"<<endl;
        flag = 0;
        return;
    }
    double r = 0;
    int row = X.getNumRow();
    for(int i = 2;i<=row;i++){
        for(int j=row;j>=i;j--){
            Y(j,1) = (Y(j,1) - Y(j-1,1))/(X(j,1) - X(j-i+1,1));
//            cout<<Y(j,1)<<endl;
        }
    }
    cout<<"差商表为:"<<endl;
    Y.printMatrix(Y);

    for(int i = 1;i<=row;i++){
        double temp = 1;
        for(int j = 1;j<i;j++){
            temp *= (a-X(j,1));
        }
        temp = temp * Y(i,1);
        r += temp;
    }
    x = r;
}