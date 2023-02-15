//
// Created by YH Q on 2022/4/1.
//
/*
 *
 * �����������ֵ������
 * ����������������
 * ���� + - * / =
 * ���� += -= *= /=
 *
 */
#ifndef CALCULATIONMETHOD_VECTOR_H
#define CALCULATIONMETHOD_VECTOR_H
#pragma once
/* ����vector.h��ͷ�ļ� */
#include "complex" //���븴��

using namespace std;

template<typename T>
class Vector {
public:
    Vector() {
        num = 0;
    }

    Vector(int n) {
        num = n;
        data = new T[n]; //����ռ�
        for (int i = 0; i < n; i++) {  //��ʼ��
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

    //��ʼ��
    void init(int n) {
        num = n;
        data = new T[n]; //����ռ�
        for (int i = 0; i < n; i++) {  //��ʼ��
            data[i] = 0.0;
        }
    }

    //�ջؿռ�
    void clear() {
        if (num > 0) {
            delete data;
            num = 0;
        }
    }


    //����һ���ӿڣ������޸����������ֵ
    int &size() {  // ���������Ĵ�С
        return num;
    }

    T &operator[](int i) {
        return data[i];
    }

    T &operator()(int i) {
        return data[i - 1];
    }

    //����+��=
    template<typename E>
    Vector<T> &operator=(const E &val) { //��ֵ��������
        for (int i = 0; i < num; i++) {
            data[i] = val;
        }
        return *this;
    }

    Vector<T> &operator=(const Vector<T> &vec) { //ͬһ�����͵�������ֵ ��Сһ��
        for (int i = 0; i < num; i++) {
            data[i] = vec.data[i];
        }
        return *this;
    }

    //˫Ŀ���������Ԫ����
    //���ؼӷ�
    friend Vector<T> operator+(const Vector<T> &vec1, const Vector<T> &vec2) {
        Vector<T> v(vec1.num);
        for (int i = 0; i <= vec1.num; ++i) {
            v.data[i] = vec1.data[i] + vec2.data[i];
        }
        return v;
    }

    //���ؼ���
    friend Vector<T> operator-(const Vector<T> &vec1, const Vector<T> &vec2) {
        Vector<T> v(vec1.num);
        for (int i = 0; i <= vec1.num; ++i) {
            v.data[i] = vec1.data[i] - vec2.data[i];
        }
        return v;
    }

    //���س˷�
    friend Vector<T> operator*(const Vector<T> &vec1, const Vector<T> &vec2) {
        Vector<T> v(vec1.num);
        for (int i = 0; i <= vec1.num; ++i) {
            v.data[i] = vec1.data[i] * vec2.data[i];
        }
        return v;
    }

    //���س���
    friend Vector<T> operator/(const Vector<T> &vec1, const Vector<T> &vec2) {
        Vector<T> v(vec1.num);
        for (int i = 0; i <= vec1.num; ++i) {
            v.data[i] = vec1.data[i] / vec2.data[i];
        }
        return v;
    }

    //����+���ֺ�����
    template<typename E>
    friend Vector<T> operator+(const E &val, const Vector<T> &vec) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] + val;
        }
        return v;
    }

    //���ؼӷ�  ����������
    template<typename E>
    friend Vector<T> operator+(const Vector<T> &vec, const E &val) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] + val;
        }
        return v;
    }

    //����-���ֺ�����
    template<typename E>
    friend Vector<T> operator-(const E &val, const Vector<T> &vec) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] - val;
        }
        return v;
    }

    //���ؼ���  ����������
    template<typename E>
    friend Vector<T> operator-(const Vector<T> &vec, const E &val) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] - val;
        }
        return v;
    }

    //����* ���ֺ�����
    template<typename E>
    friend Vector<T> operator*(const E &val, const Vector<T> &vec) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] * val;
        }
        return v;
    }

    //���ؼ���  ����������
    template<typename E>
    friend Vector<T> operator*(const Vector<T> &vec, const E &val) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] * val;
        }
        return v;
    }


    //���� /  ����/����
    template<typename E>
    friend Vector<T> operator/(const Vector<T> &vec, const E &val) {
        Vector<T> v(vec.num);
        for (int i = 0; i < vec.num; ++i) {
            v.data[i] = vec.data[i] / val;
        }
        return v;
    }

    //����+=
    Vector<T> &operator+=(const Vector<T> &vec) { //ͬһ�����͵�������ֵ ��Сһ��
        for (int i = 0; i < num; i++) {
            data[i] += vec.data[i];
        }
        return *this;
    }

    //����-=
    Vector<T> &operator-=(const Vector<T> &vec) { //ͬһ�����͵�������ֵ ��Сһ��
        for (int i = 0; i < num; i++) {
            data[i] -= vec.data[i];
        }
        return *this;
    }

    //����*=
    Vector<T> &operator*=(const Vector<T> &vec) { //ͬһ�����͵�������ֵ ��Сһ��
        for (int i = 0; i < num; i++) {
            data[i] *= vec.data[i];
        }
        return *this;
    }

    //���� /=
    Vector<T> &operator/=(const Vector<T> &vec) { //ͬһ�����͵�������ֵ ��Сһ��
        for (int i = 0; i < num; i++) {
            data[i] /= vec.data[i];
        }
        return *this;
    }

    //���� += ����ֵ֮������
    template<typename E>
    Vector<T> &operator+=(const E &val) { //��ֵ��������
        for (int i = 0; i < num; i++) {
            data[i] += val;
        }
        return *this;
    }

    //���� -= ����ֵ����
    template<typename E>
    Vector<T> &operator-=(const E &val) { //��ֵ��������
        for (int i = 0; i < num; i++) {
            data[i] -= val;
        }
        return *this;
    }

    //���� *= ����ֵ����
    template<typename E>
    Vector<T> &operator*=(const E &val) { //��ֵ��������
        for (int i = 0; i < num; i++) {
            data[i] *= val;
        }
        return *this;
    }

    //���� /= ����ֵ����
    template<typename E>
    Vector<T> &operator/=(const E &val) { //��ֵ��������
        for (int i = 0; i < num; i++) {
            data[i] /= val;
        }
        return *this;
    }

protected:
    int num; //�������Ⱥʹ�С
    T *data; //����������
};



#endif //CALCULATIONMETHOD_VECTOR_H
