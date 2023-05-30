/*
ЗАДАЧА 12:
Дан базис 
e1 = i;
e2 = 2j;
e3 = j + k;
Даны дважды контравариантные компоненты тензоров
P = 
4 3 -1
0 2 2
-3 4 5

T =
6 2 7
-7 5 -1
0 4 3

Найти полное и двойное скалярное произведение, проверить свойства
*/

#include <iostream>
#include "classes.cpp"

using namespace std;

int main() {
    // Дано
    Vector e1 = Vector(1, 0, 0);
    Vector e2 = Vector(0, 2, 0);
    Vector e3 = Vector(0, 1, 1);
    Basis e = Basis(e1, e2, e3);
    double **Pcomps = new double*[3];
    Pcomps[0] = new double[3] {4, 3, -1};
    Pcomps[1] = new double[3] {0, 2, 2};
    Pcomps[2] = new double[3] {-3, 4, 5};
    Tensor_II P = Tensor_II(e, 0, Pcomps);
    double **Tcomps = new double*[3];
    Tcomps[0] = new double[3] {6, 2, 7};
    Tcomps[1] = new double[3] {-7, 5, -1};
    Tcomps[2] = new double[3] {0, 4, 3};
    Tensor_II T = Tensor_II(e, 0, Tcomps);
    // Решение
    cout << P.double_scalar_product(T) << endl;
    cout << T.double_scalar_product(P) << endl;
    cout << P.total_scalar_product(T.transpose()) << endl;
    cout << P.transpose().total_scalar_product(T) << endl;
    cout << P.total_scalar_product(T) << endl;
    return 0;
}