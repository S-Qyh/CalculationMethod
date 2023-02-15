//
// Created by YH Q on 2022/4/1.
//
#pragma once
#ifndef CALCULATIONMETHOD_MATRIX_H
#define CALCULATIONMETHOD_MATRIX_H

//#include "Vector.h"
#include "iostream"
#include "math.h"
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
        assert(this->data!=NULL);
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

    //矩阵的转置。后置运算 ～

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
    double colVectorNorm_3(Matrix<T> &mat) {
        double max = fabs(mat(1, 1));
        for (int i = 1; i <= mat.num_Row; i++) {
            if (fabs(mat(i, 1)) > max) {
                max = fabs(mat(i, 1));
            }
        }
        return max;
    }

    //2 ?????  ????A????ù???????????A??????????????????????



    void printMatrix(Matrix<T> &mat) const {
        for (int i = 1; i <= mat.num_Row; i++) {
            for (int j = 1; j <= mat.num_Col; j++) {
                cout << mat(i, j) << " ";
            }
            cout << endl;
        }
        cout << endl;
    }


    void inputMat(Matrix<T> &mat) {
        for (int i = 1; i <= this->num_Row; i++) {
            for (int j = 1; j <= this->num_Col; j++) {
                cin >> mat(i, j);
            }
        }
    }

    //矩阵范数还没实现
    void MatrixNorm(int type,int &flag,double &result){

    }


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
