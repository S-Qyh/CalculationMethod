//
// Created by YH Q on 2022/4/1.
//
#pragma once
#ifndef CALCULATIONMETHOD_MATRIX_H
#define CALCULATIONMETHOD_MATRIX_H

#include "Vector.h"
#include "iostream"
using namespace std;

template<typename T>
class Matrix {
public:
    Matrix() {
        num_Row = 0;
        num_Col = 0;
    }

    //???????n*n?????
    Matrix(int n) {
        num_Row = n;
        num_Col = n;
        data = new T[n * n];
        for (int i = 0; i < n * n; ++i) {
            data[i] = 0.0;
        }
    }

    //?????к??У???????????????к??ж???0?????
    Matrix(int m, int n) {
        num_Row = m;
        num_Col = n;
        data = new T[m * n];
        for (int i = 0; i < m * n; ++i) {
            data[i] = 0.0;
        }
    }

    //????????
    Matrix(Matrix<T> &mat) {
        num_Row = mat.num_Row;
        num_Col = mat.num_Col;
        data = new T[num_Row * num_Col];
        for (int i = 0; i < num_Row * num_Col; ++i) {
            data[i] = mat.data[i];
        }
    }

    ~Matrix() {
        if (num_Row * num_Col > 0) {
            delete data;
            num_Row = 0;
            num_Col = 0;
        }
    }

    //??????
    void init(int m, int n) {
        num_Row = m;
        num_Col = n;
        data = new T[m * n];
        for (int i = 0; i < m * n; ++i) {
            data[i] = 0.0;
        }
    }

    //?????
    void clear() {
        if (num_Row * num_Col > 0) {
            delete data;
            num_Row = 0;
            num_Col = 0;
        }
    }

    //????init??????????????? ????n*n?????
    void inti(int n) {
        num_Row = n;
        num_Col = n;
        data = new T[n * n];
        for (int i = 0; i < n * n; ++i) {
            data[i] = 0.0;
        }
    }

    //???????
    int &nRow() {
        return num_Row; //??????????
    }

    int &nCol() {
        return num_Col;//??????????
    }

    //?????????????? ????????
    T &operator[](int i) {
        return data[i];
    }

    //?????i?е?j?????
    T &operator()(int i, int j) {
        //??????   ???????n??  ??j - 1?? * num_Row ??????ж???????????????е?λ??
        return data[i - 1 + (j - 1) * num_Row];
    }

    //?????????
    //???? = ???????????
    template<typename E>
    Matrix<T> &operator=(const E &val) { //???????????
        for (int i = 0; i < num_Row * num_Col; i++) {
            data[i] = val;
        }
        return *this;
    }

    //??????????
    Matrix<T> &operator=(const Matrix<T> &mat) { //????????????? ??С???
        for (int i = 0; i < num_Row * num_Col; i++) {
            data[i] = mat.data[i];
        }
        return *this;
    }

    //???? +=
    //????????
    template<typename E>
    Matrix<T> &operator+=(const E &val) {
        for (int i = 0; i < num_Row * num_Col; i++) {
            data[i] += val;
        }
        return *this;
    }

    //???????? +=  ????????????
    Matrix<T> &operator+=(const Matrix<T> &mat) {
        for (int i = 0; i < num_Row * num_Col; i++) {
            data[i] += mat.data[i];
        }
        return *this;
    }

    // ???? -= ?????????
    template<typename E>
    Matrix<T> &operator-=(const E &val) {
        for (int i = 0; i < num_Row * num_Col; i++) {
            data[i] -= val;
        }
        return *this;
    }

    // -= ?????????
    Matrix<T> &operator-=(const Matrix<T> &mat) {
        for (int i = 0; i < num_Row * num_Col; i++) {
            data[i] -= mat.data[i];
        }
        return *this;
    }

    //????*= ?????????
    template<typename E>
    Matrix<T> &operator*=(const E &val) {
        for (int i = 0; i < num_Row * num_Col; i++) {
            data[i] *= val;
        }
        return *this;
    }

    //????*= ?????????
    Matrix<T> &operator*=(const Matrix<T> &mat) {
        for (int i = 0; i < num_Row * num_Col; i++) {
            data[i] *= mat.data[i];
        }
        return *this;
    }

    //???? /= ?????????
    template<typename E>
    Matrix<T> &operator/=(const E &val) {
        for (int i = 0; i < num_Row * num_Col; i++) {
            data[i] /= val;
        }
        return *this;
    }

    //???? /= ?????????
    Matrix<T> &operator/=(const Matrix<T> &mat) {
        for (int i = 0; i < num_Row * num_Col; i++) {
            data[i] /= mat.data[i];
        }
        return *this;
    }

    //???? + ????????
    template<typename E>
    friend Matrix<T> operator+(const E &val, const Matrix<T> &mat) {
        Matrix<T> v(mat.num_Row, mat.num_Col);
        for (int i = 0; i < mat.num_Row * mat.num_Col; i++) {
            v.data[i] = val + mat.data[i];
        }
        return v;
    }

    // ???? + ???
    template<typename E>
    friend Matrix<T> operator+(const Matrix<T> &mat, const E &val) {
        Matrix<T> v(mat.num_Row, mat.num_Col);
        for (int i = 0; i < mat.num_Row * mat.num_Col; i++) {
            v.data[i] = mat.data[i] + val;
        }
        return v;
    }

    //???? + ????   ???????????
    friend Matrix<T> operator+(const Matrix<T> &mat1, const Matrix<T> &mat2) {
        Matrix<T> v(mat1.num_Row, mat1.num_Col);
        for (int i = 0; i < mat1.num_Row * mat1.num_Col; ++i) {
            v.data[i] = mat1.data[i] + mat2.data[i];
        }
        return v;
    }

    //???? -    ??? - ????
    template<typename E>
    friend Matrix<T> operator-(const E &val, const Matrix<T> &mat) {
        Matrix<T> v(mat.num_Row, mat.num_Col);
        for (int i = 0; i < mat.num_Row * mat.num_Col; i++) {
            v.data[i] = val - mat.data[i];
        }
        return v;
    }

    //???? -  ???? - ???

    template<typename E>
    friend Matrix<T> operator-(const Matrix<T> &mat, const E &val) {
        Matrix<T> v(mat.num_Row, mat.num_Col);
        for (int i = 0; i < mat.num_Row * mat.num_Col; i++) {
            v.data[i] = mat.data[i] - val;
        }
        return v;
    }

    //???? - ????
    friend Matrix<T> operator-(const Matrix<T> &mat1, const Matrix<T> &mat2) {
        Matrix<T> v(mat1.num_Row, mat1.num_Col);
        for (int i = 0; i < mat1.num_Row * mat1.num_Col; ++i) {
            v.data[i] = mat1.data[i] - mat2.data[i];
        }
        return v;
    }

    //??? * ????
    template<typename E>
    friend Matrix<T> operator*(const E &val, const Matrix<T> &mat) {
        Matrix<T> v(mat.num_Row, mat.num_Col);
        for (int i = 0; i < mat.num_Row * mat.num_Col; ++i) {
            v.data[i] = val * mat.data[i];
        }
        return v;
    }

    //???? *  ???
    template<typename E>
    friend Matrix<T> operator*(const Matrix<T> &mat, const E &val) {
        Matrix<T> v(mat.num_Row, mat.num_Col);
        for (int i = 0; i < mat.num_Row * mat.num_Col; i++) {
            v.data[i] = mat.data[i] * val;
        }
        return v;
    }

    //???? * ????
    friend Matrix<T> operator*(const Matrix<T> &mat1, const Matrix<T> &mat2) {
        Matrix<T> v(mat1.num_Row, mat1.num_Col);
        for (int i = 0; i < mat1.num_Row * mat1.num_Col; i++) {
            v.data[i] = mat1.data[i] * mat2.data[i];
        }
        return v;
    }

    //???? / ????
    friend Matrix<T> operator/(const Matrix<T> &mat1, const Matrix<T> &mat2) {
        Matrix<T> v(mat1.num_Row, mat1.num_Col);
        for (int i = 0; i < mat1.num_Row * mat1.num_Col; ++i) {
            v.data[i] = mat1.data[i] / mat2.data[i];
        }
        return v;
    }

    //  ???? / ????
    template<typename E>
    friend Matrix<T> operator/(const Matrix<T> &mat, const E &val) {
        Matrix<T> v(mat.num_Row, mat.num_Col);
        for (int i = 0; i < mat.num_Row * mat.num_Col; ++i) {
            v.data[i] = mat.data[i] / val;
        }
        return v;
    }

    //1-???? ?к???????????о??????????????????????
    double vectorNorm_1(Matrix<T> &mat) {
        int row = mat.num_Row;
        int col = mat.num_Col;
        double sum = 0;
        double temp[col];
        for (int i = 1; i <= col; i++) {
            for (int j = 1; j <= row; j++) {
                sum = sum + fabs(mat(j, i));
            }
            temp[i] = sum;
            sum = 0;
        }

        for (int i = 0; i < col; ++i) {
            if (temp[0] > temp[i]) {
                continue;
            } else {
                temp[0] = temp[i + 1];
            }
        }
        return temp[0];
    }

    //??????
    double vectorNorm_3(Matrix<T> &mat) {
        double max = fabs(mat(1, 1));
        for (int i = 1; i <= mat.num_Row; i++) {
            if (fabs(mat(i, 1)) > max) {
                max = fabs(mat(i, 1));
            }
        }
        return max;
    }

    //2 ?????  ????A????ù???????????A??????????????????????


    //????????
//    void Gauss(Matrix<T> &mat, Matrix<T> &b, int n) {
//        for (int i = 1; i <= n; i++) {
//            int colnum = i;
//            double max = 0;
//            for (int k = i; k < n; k++) {
//                if (fabs(mat(k, i) > max)) {
//                    max = fabs(mat(k, i));
//                    colnum = i;
//                }
//            }
//
//            for (int w = i; w < n; w++) {
//                double temp = mat(i, w);
//                mat(i, w) = mat(colnum, w);
//                mat(colnum, w) = temp;
//            }
//
//            cout << "??" << i << "???????\n";
//            mat.printMatrix(mat);
//
//
//            for (int j = i + 1; j <= n; j++) {
//                double M = mat(j, i) / mat(i, i);
////                cout << M << endl;
//                for (int k = 1; k < n + 1; k++) {
//                    mat(j, k) = mat(j, k) - M * mat(i, k);
////                    cout << mat(j, k) << "  ";
//                }
//                b(j, 1) = b(j, 1) - M * b(i, 1);
////                cout << b(j, 1) << endl;
//            }
//        }
//        cout << "????????????\n";
//        mat.printMatrix(mat);
//
//        double x[n];
//        x[n] = b(n, 1) / mat(n, n);
//        cout << "x" << n << "???" << x[n] << endl;
//        for (int i = n - 1; i >= 1; i--) {
//
//            double sum = 0;
//            for (int k = i + 1; k <= n; k++) {
//                sum = sum + mat(i, k) * x[k];
//            }
//            x[i] = (b(i, 1) - sum) / mat(i, i);
//            cout << "x" << i << "???" << x[i] << endl;
//        }
//    }

    void printMatrix(Matrix<T> &mat) const {
        for (int i = 1; i <= mat.num_Row; i++) {
            for (int j = 1; j <= mat.num_Col; j++) {
                cout << mat(i, j) << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

//    void matrixLU(Matrix<T> &mat, Matrix<T> &b, int n) {
//        Matrix<T> L(mat.num_Row, mat.num_Col);
//        Matrix<T> U(mat.num_Row, mat.num_Col);
//
//        for (int i = 1; i <= n; i++) {
//            U(1, i) = mat(1, i);
//        }
//
//        for (int i = 1; i <= n; i++) {
//            L(i, i) = 1;
//            L(i, 1) = mat(i, 1) / U(1, 1);
//        }
////        cout<<U(1,1)<<U(1,2)<<U(1,3)<<endl;
////        cout<<L(1,1)<<L(2,1)<<L(3,1)<<endl;
//        for (int t = 1; t <= n; t++) {
//            for (int i = t; i < n + 1; i++) {
//                double sum_LU = 0;
//                for (int k = 1; k < t; k++) {
//                    sum_LU = sum_LU + L(t, k) * U(k, i);
//                }
//                U(t, i) = mat(t, i) - sum_LU;
//                sum_LU = 0;
//                for (int k = 1; k < t; k++) {
//                    sum_LU = sum_LU + L(i, k) * U(k, t);
//                }
//                L(i, t) = (mat(i, t) - sum_LU) / U(t, t);
//            }
//        }
//        cout << "L???????\n";
//        L.printMatrix(L);
//
//        cout << "U???????\n";
//        U.printMatrix(U);
//
//        double y[n];
//        y[1] = b(1, 1);
//        for (int i = 2; i <= n; i++) {
//            double sum = 0;
//            for (int k = 1; k <= i - 1; k++) {
//                sum = sum + L(i, k) * y[k];
//            }
//            y[i] = b(i, 1) - sum;
//        }
//        for (int i = 1; i <= b.num_Row; i++)
//            cout << "y" << i << "???" << y[i] << endl;
//        cout << endl;
//
//        double x[n];
//        x[n] = y[n] / U(n, n);
//
//        for (int i = n - 1; i >= 1; i--) {
//            double sum = 0;
//            for (int k = i + 1; k <= n; k++) {
//                sum = sum + U(i, k) * x[k];
//            }
//            x[i] = (y[i] - sum) / U(i, i);
//        }
//        for (int i = 1; i <= b.num_Row; i++)
//            cout << "x" << i << "???" << x[i] << endl;
//    }

    void inputMat(Matrix<T> &mat) {
        for (int i = 1; i <= this->num_Row; i++) {
            for (int j = 1; j <= this->num_Col; j++) {
                cin >> mat(i, j);
            }
        }
    }

//    void jacobiIterator(Matrix<T> &A, Matrix<T> &B, double e, int M, int flag, Matrix<T> X0, Matrix<T> &X) {
//        int n = A.num_Row;
//        Matrix<T> temp(n,1);
//        double sum = 0;
//        int k = 1;
//        while (true) {
//            for (int i = 1; i <= n; i++) {
//                for (int j = 1; j <= n; j++) {
//                    if(j == i){
//                        continue;
//                    }
//                    sum += A(i,j) * X0(j,1);
//                }
//
//                X(i,1) = (B(i,1) - sum) / A(i,i);
//
//                sum = 0;
//            }
//            cout<<"??"<<k<<"?ν????\n";
//            X.printMatrix(X);
//            temp = X - X0;
////            cout<<X.vectorNorm_3(temp)<<" ";
//            if(X.vectorNorm_3(temp) < e){
//                cout<<endl;
//                cout<<"???????\n";
//                X.printMatrix(X);
//                break;
//            } else{
//                X0 = X;
//            }
//            if(k == M){
//                cout<<M<<"?ε???δ??С??e??????";
//                break;
//            } else{
//                k++;
//            }
//
//        }
//
//    }
    int getNumRow() const {
        return num_Row;
    }

    int getNumCol() const {
        return num_Col;
    }


protected:
    int num_Row;//??
    int num_Col;//??
    T *data; //?б???
};


#endif //CALCULATIONMETHOD_MATRIX_H
