//#include "Vector.h"
#include "iostream"
#include "Matrix.h"
#include "LinearEquationsSolve.h"

using namespace std;

int main(int argc, const char *argv[]) {
    Matrix<double> a(4, 4), b(4, 4), c(3, 3),x(4,1);
    LinearEquationsSolve le;
    a(1,1) = 2;
    a(1,2) = -1;
    a(1,3) = 3;
    a(1,4) = 2,

    a(2,1 ) = 4;
    a(2,2) = 2;
    a(2,3) = 5;
    a(2,4) = 5;

    a(3,1) = 1;
    a(3,2) = 2;
    a(3,3) = 0;
    a(3,4) = 8;

    a(4,1) = 1;
    a(4,2) = -3;
    a(4,3) = 5;
    a(4,4) = 2;

    b(1,1) = 1;
    b(2,1) = 4;
    b(3,1) = 7;
    b(4,1) = 3;

    int flag = 0;

    le.Gauss(a,b,x,flag);

    le.GaussCol(a,b,flag,x);

}