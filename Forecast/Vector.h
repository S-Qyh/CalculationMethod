//
// Created by YH Q on 2022/4/1.
//
/*
 *
 * 完成向量和数值的运算
 * 向量和向量的运算
 * 重载 + - * / =
 * 重载 += -= *= /=
 *
 */
#ifndef CALCULATIONMETHOD_VECTOR_H
#define CALCULATIONMETHOD_VECTOR_H
#pragma once
/* 创建vector.h的头文件 */
#include "complex" //加入复数

using namespace std;

template<typename T>
class Vector {
public:
    Vector() {
        num = 0;
    }

    Vector(int n) {
        num = n;
        data = new T[n]; //分配空间
        for (int i = 0; i < n; i++) {  //初始化
            data[i] = 0.0;
        }
    }

    Vector(const Vector<T> &vec) {
        num = vec.num;
        data = new T[num];
        for (int i = 0; i < num; i++) {
            data[i] = vec.data[i];
        }
    }

    ~Vector() {
        if (num > 0) {
            delete data;
            num = 0;
        }
    }

    //初始化
    void init(int n) {
        num = n;
        data = new T[n]; //分配空间
        for (int i = 0; i < n; i++) {  //初始化
            data[i] = 0.0;
        }
    }

    //收回空间
    void clear() {
        if (num > 0) {
            delete data;
            num = 0;
        }
    }


    //创建一个接口，用来修改向量里面的值
    int &size() {  // 返回向量的大小
        return num;
    }

    T &operator[](int i) {
        return data[i];
    }

    T &operator()(int i) {
        return data[i - 1];
    }

    //重载+和=
    template<typename E>
    Vector<T> &operator=(const E &val) { //数值赋给向量
        for (int i = 0; i < num; i++) {
            data[i] = val;
        }
        return *this;
    }

    Vector<T> &operator=(const Vector<T> &vec) { //同一种类型的向量赋值 大小一样
        for (int i = 0; i < num; i++) {
            data[i] = vec.data[i];
        }
        return *this;
    }

    //双目运算符用友元方便
    //重载加法
    friend Vector<T> operator+(const Vector<T> &vec1, const Vector<T> &vec2) {
        Vector<T> v(vec1.num);
        for (int i = 0; i <= vec1.num; ++i) {
            v.data[i] = vec1.data[i] + vec2.data[i];
        }
        return v;
    }

    //重载减法
    friend Vector<T> operator-(const Vector<T> &vec1, const Vector<T> &vec2) {
        Vector<T> v(vec1.num);
        for (int i = 0; i <= vec1.num; ++i) {
            v.data[i] = vec1.data[i] - vec2.data[i];
        }
        return v;
    }

    //重载乘法
    friend Vector<T> operator*(const Vector<T> &vec1, const Vector<T> &vec2) {
        Vector<T> v(vec1.num);
        for (int i = 0; i <= vec1.num; ++i) {
            v.data[i] = vec1.data[i] * vec2.data[i];
        }
        return v;
    }

    //重载除法
    friend Vector<T> operator/(const Vector<T> &vec1, const Vector<T> &vec2) {
        Vector<T> v(vec1.num);
        for (int i = 0; i <= vec1.num; ++i) {
            v.data[i] = vec1.data[i] / vec2.data[i];
        }
        return v;
    }

    //重载+数字和向量
    template<typename E>
    friend Vector<T> operator+(const E &val, const Vector<T> &vec) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] + val;
        }
        return v;
    }

    //重载加法  向量和数字
    template<typename E>
    friend Vector<T> operator+(const Vector<T> &vec, const E &val) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] + val;
        }
        return v;
    }

    //重载-数字和向量
    template<typename E>
    friend Vector<T> operator-(const E &val, const Vector<T> &vec) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] - val;
        }
        return v;
    }

    //重载减法  向量和数字
    template<typename E>
    friend Vector<T> operator-(const Vector<T> &vec, const E &val) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] - val;
        }
        return v;
    }

    //重载* 数字和向量
    template<typename E>
    friend Vector<T> operator*(const E &val, const Vector<T> &vec) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] * val;
        }
        return v;
    }

    //重载减法  向量和数字
    template<typename E>
    friend Vector<T> operator*(const Vector<T> &vec, const E &val) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] * val;
        }
        return v;
    }


    //重载 /  向量/数字
    template<typename E>
    friend Vector<T> operator/(const Vector<T> &vec, const E &val) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] / val;
        }
        return v;
    }

    //重载+=
    Vector<T> &operator+=(const Vector<T> &vec) { //同一种类型的向量赋值 大小一样
        for (int i = 0; i < num; i++) {
            data[i] += vec.data[i];
        }
        return *this;
    }

    //重载-=
    Vector<T> &operator-=(const Vector<T> &vec) { //同一种类型的向量赋值 大小一样
        for (int i = 0; i < num; i++) {
            data[i] -= vec.data[i];
        }
        return *this;
    }

    //重载*=
    Vector<T> &operator*=(const Vector<T> &vec) { //同一种类型的向量赋值 大小一样
        for (int i = 0; i < num; i++) {
            data[i] *= vec.data[i];
        }
        return *this;
    }

    //重载 /=
    Vector<T> &operator/=(const Vector<T> &vec) { //同一种类型的向量赋值 大小一样
        for (int i = 0; i < num; i++) {
            data[i] /= vec.data[i];
        }
        return *this;
    }

    //重载 += 和数值之间运算
    template<typename E>
    Vector<T> &operator+=(const E &val) { //数值赋给向量
        for (int i = 0; i < num; i++) {
            data[i] += val;
        }
        return *this;
    }

    //重载 -= 和数值运算
    template<typename E>
    Vector<T> &operator-=(const E &val) { //数值赋给向量
        for (int i = 0; i < num; i++) {
            data[i] -= val;
        }
        return *this;
    }

    //重载 *= 和数值运算
    template<typename E>
    Vector<T> &operator*=(const E &val) { //数值赋给向量
        for (int i = 0; i < num; i++) {
            data[i] *= val;
        }
        return *this;
    }

    //重载 /= 和数值运算
    template<typename E>
    Vector<T> &operator/=(const E &val) { //数值赋给向量
        for (int i = 0; i < num; i++) {
            data[i] /= val;
        }
        return *this;
    }

protected:
    int num; //向量长度和大小
    T *data; //向量的数据
};



#endif //CALCULATIONMETHOD_VECTOR_H
