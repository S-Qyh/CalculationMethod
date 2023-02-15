//
// Created by YH Q on 2022/4/11.
//

#include "LinearEquationsSolve.h"
#include <cmath>
#include <iostream>

/*
 * 修改内容：
 *  1 函数通过参数回代结果
 *  2 将LU分解的过程单独拿出来
 *  3 列主元消去法判断是否需要  换行  加一条  if  语句
 *
 */


using namespace std;

void LinearEquationsSolve::GaussCol(Matrix<double> &A, Matrix<double> &B,int &flag,Matrix<double> &x) {
    int n = A.getNumCol();
    for (int i = 1; i <= n ; i++) {
        ColPivot(A,B,i);
        int colnum = i;
        double max = 0;
        for (int k = i; k < n; k++) {
            if (fabs(A(i,i))<0.00001){
                flag = 0;
                return;
            }
            if (fabs(A(k, i)) > max) {
                max = fabs(A(k, i));
                colnum = i;
            }
        }

        cout << "第" << i << "次消元：\n";
        A.printMatrix(A);


        for (int j = i + 1; j <= n; j++) {
            double M = A(j, i) / A(i, i);
//                cout << M << endl;
            for (int k = 1; k < n + 1; k++) {
                A(j, k) = A(j, k) - M * A(i, k);
//                    cout << A(j, k) << "  ";
            }
            B(j, 1) = B(j, 1) - M * B(i, 1);
//                cout << B(j, 1) << endl;
        }
    }
    cout << "高斯消元后的矩阵：\n";
    A.printMatrix(A);

    x(n,1) = B(n, 1) / A(n, n);
    cout << "x" << n << "解为" << x(n,1) << endl;
    for (int i = n - 1; i >= 1; i--) {

        double sum = 0;
        for (int k = i + 1; k <= n; k++) {
            sum = sum + A(i, k) * x(k,1);
        }
        x(i,1) = (B(i, 1) - sum) / A(i, i);
        cout << "x" << i << "解为" << x(i,1) << endl;
        flag = 1;
    }
}

LinearEquationsSolve::LinearEquationsSolve() {

}

void LinearEquationsSolve::matrixLU(Matrix<double> &mat, Matrix<double> &b, int n) {
    Matrix<double> L(mat.getNumRow(), mat.getNumCol());
    Matrix<double> U(mat.getNumRow(), mat.getNumCol());

    for (int i = 1; i <= n; i++) {
        U(1, i) = mat(1, i);
    }

    for (int i = 1; i <= n; i++) {
        L(i, i) = 1;
        L(i, 1) = mat(i, 1) / U(1, 1);
    }
//        cout<<U(1,1)<<U(1,2)<<U(1,3)<<endl;
//        cout<<L(1,1)<<L(2,1)<<L(3,1)<<endl;
    for (int t = 1; t <= n; t++) {
        for (int i = t; i < n + 1; i++) {
            double sum_LU = 0;
            for (int k = 1; k < t; k++) {
                sum_LU = sum_LU + L(t, k) * U(k, i);
            }
            U(t, i) = mat(t, i) - sum_LU;
            sum_LU = 0;
            for (int k = 1; k < t; k++) {
                sum_LU = sum_LU + L(i, k) * U(k, t);
            }
            L(i, t) = (mat(i, t) - sum_LU) / U(t, t);
        }
    }
    cout << "L矩阵为：\n";
    L.printMatrix(L);

    cout << "U矩阵为：\n";
    U.printMatrix(U);

    double y[n];
    y[1] = b(1, 1);
    for (int i = 2; i <= n; i++) {
        double sum = 0;
        for (int k = 1; k <= i - 1; k++) {
            sum = sum + L(i, k) * y[k];
        }
        y[i] = b(i, 1) - sum;
    }
    for (int i = 1; i <= b.getNumRow(); i++)
        cout << "y" << i << "解为" << y[i] << endl;
    cout << endl;

    double x[n];
    x[n] = y[n] / U(n, n);

    for (int i = n - 1; i >= 1; i--) {
        double sum = 0;
        for (int k = i + 1; k <= n; k++) {
            sum = sum + U(i, k) * x[k];
        }
        x[i] = (y[i] - sum) / U(i, i);
    }
    for (int i = 1; i <= b.getNumRow(); i++)
        cout << "x" << i << "解为" << x[i] << endl;
}

void
LinearEquationsSolve::jacobiIterator(Matrix<double> &A, Matrix<double> &B, double e, int M, int &flag,
                                     Matrix<double> X0,
                                     Matrix<double> &X) {
    int n = A.getNumRow();
    Matrix<double> temp(n, 1);
    double sum = 0;
    int k = 1;
    while (true) {
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                if (j == i) {
                    continue;
                }
                sum += A(i, j) * X0(j, 1);
            }

            X(i, 1) = (B(i, 1) - sum) / A(i, i);

            sum = 0;
        }
        cout << "第" << k << "次结果：\n";
        X.printMatrix(X);
        temp = X - X0;
//            cout<<X.colVectorNorm_3(temp)<<" ";
        if (X.colVectorNorm_3(temp) < e) {
            cout << endl;
            cout << "最后结果：\n";
            X.printMatrix(X);
            flag = 1;
            break;
        } else {
            X0 = X;
        }
        if (k == M) {
            cout << M << "次迭代未达到小于e的结果。";
            flag = 0;
            break;
        } else {
            k++;
        }

    }

}

void LinearEquationsSolve::gaussSD(Matrix<double> &A, Matrix<double> &B, double e, int M, int &flag, Matrix<double> X0,
                                   Matrix<double> &X) {
    int n = A.getNumRow();
    Matrix<double> temp(n, 1);
    Matrix<double> temp1(n, 1);
    temp1 = X0;
    int k = 1;
    while (true) {
        for (int i = 1; i <= n; i++) {
            double sum = 0;
            for (int j = 1; j <= n; j++) {
                if (j == i) {
                    continue;
                }
                sum += A(i, j) * X0(j, 1);
            }
            X(i, 1) = (B(i, 1) - sum) / A(i, i);
            X0(i, 1) = X(i, 1);
//            sum = 0;
        }
        cout << "第" << k << "次结果：\n";
        X.printMatrix(X);
        temp = X - temp1;
        temp1 = X0;
//            cout<<X.colVectorNorm_3(temp)<<" ";
        if (X.colVectorNorm_3(temp) < e) {
            cout << endl;
            cout << "最后结果：\n";
            X.printMatrix(X);
            flag = 1;
            break;
        } else {
            X0 = X;
        }
        if (k == M) {
            cout << M << "次迭代未达到小于e的结果。";
            flag = 0;
            break;
        } else {
            k++;
        }

    }

}

void LinearEquationsSolve::Gauss(Matrix<double> &A, Matrix<double> &B,Matrix<double> &x,int &flag) {
    int n = A.getNumCol();
    for (int i = 1; i <= n ; i++) {
        int colnum = i;
        double max = 0;
        for (int k = i; k < n; k++) {
            if (fabs(A(i,i))<0.00001){
                flag = 0;
                return;
            }
            if (fabs(A(k, i)) > max) {
                max = fabs(A(k, i));
                colnum = i;
            }
        }

        cout << "第" << i << "次消元：\n";
        A.printMatrix(A);


        for (int j = i + 1; j <= n; j++) {
            double M = A(j, i) / A(i, i);
//                cout << M << endl;
            for (int k = 1; k < n + 1; k++) {
                A(j, k) = A(j, k) - M * A(i, k);
//                    cout << A(j, k) << "  ";
            }
            B(j, 1) = B(j, 1) - M * B(i, 1);
//                cout << B(j, 1) << endl;
        }
    }
    cout << "高斯消元后的矩阵：\n";
    A.printMatrix(A);

    x(n,1) = B(n, 1) / A(n, n);
    cout << "x" << n << "解为" << x(n,1) << endl;
    for (int i = n - 1; i >= 1; i--) {

        double sum = 0;
        for (int k = i + 1; k <= n; k++) {
            sum = sum + A(i, k) * x(k,1);
        }
        x(i,1) = (B(i, 1) - sum) / A(i, i);
        cout << "x" << i << "解为" << x(i,1) << endl;
        flag = 1;
    }
}

void LinearEquationsSolve::ColPivot(Matrix<double> &A, Matrix<double> &B,int k) {
    int row = A.getNumRow();
    int col = A.getNumCol();

    double max = fabs(A(k,k));
    int pos = k;
    for (int i = k + 1;i <= row;i++){
        if(max< fabs(A(i,k))){
            pos = i;
            max = fabs(A(i,k));
        }
    }
    if(pos>k){
        for(int j = 1;j<=col;j++){
            double temp = A(k,j);
            A(k,j) = A(pos,j);
            A(pos,j) = temp;
        }
        double temp = B(k,1);
        B(k,1) = B(pos,1);
        B(pos,1) = temp;
    }
}

