#include <iostream>

#include "Matrix.hpp"

double idenity(double x) { return x; }

int main() {
    Matrix a;           //Empty matrix, no data (maybe can resize in futere verions?)
    Matrix b(3, 3);     //A 3x3 Matrix with 0 as a inital value
    Matrix c(3, 3, 1);  // A 3x3 Matrix with 1 as initial value
    Matrix d(3, 4, 5);  // A 4x3 Matrix with 5 as inital value

    b *= c;         // Mutiply two Matrixes (mutiple matrix not suported yet ex: b *= c + c)
    c += b;         //  Sum two Matrixes (mutiple matrix not suported yet ex: c += c + c)
    b = c + b - c;  // Atribuition works (not mix +/- with *);
    b = c * d;      // For now i dont recomend make mutiple operations with mutiplication in the same statement

    //Scalar operations

    b = b + 1;
    b = b - 5;
    b = b * 10;
    b = b / 2;

    b += 3;
    b *= 0.5;
    b -= 4;
    b /= 2;

    b.map(idenity); // Map apply a function in each element of the Matrix

    a.print();  // is equivalent a std::cout << a;

    std::cout
        << "Matrix A:\n"
        << a
        << "Matrix B:\n"
        << b
        << "Matrix C:\n"
        << c
        << "Matrix D:\n"
        << d
        << std::endl;

    std::cout << "[PASS]" << std::endl;
    return 0;
}
