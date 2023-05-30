/*
Дан ОНБ базис 
e1 = i;
e2 = j;
e3 = k;
Даны компоненты T - тензора 2 ранга:
1 7 3
2 9 0
3 -4 4
Найти его сопутствующий вектор.
Даны компоненты A - антисимметричного тензора 2 ранга:
0 8 -4
-8 0 6
4 -6 0 
Найти его аксиальный вектор
*/

#include <iostream>
#include "classes.cpp"

using namespace std;

int main() {
    // Дано
    Vector e1 = Vector(1, 0, 0);
    Vector e2 = Vector(0, 1, 0);
    Vector e3 = Vector(0, 0, 1);
    Basis e = Basis(e1, e2, e3);
    double **Tcomps = new double*[3];
    Tcomps[0] = new double[3] {1, 7, 3};
    Tcomps[1] = new double[3] {2, 9, 0};
    Tcomps[2] = new double[3] {3, -4, 4};
    Tensor_II T = Tensor_II(e, 0, Tcomps);
    double **Acomps = new double*[3];
    Acomps[0] = new double[3] {0, 8, -4};
    Acomps[1] = new double[3] {-8, 0, 6};
    Acomps[2] = new double[3] {4, -6, 0};
    Tensor_II A = Tensor_II(e, 0, Acomps);
    A.show(0);
    // Решение
    Tensor_I t = T.double_scalar_product(LeviCivita()) * 0.5;
    cout << "t = " << endl;
    t.show(true);
    Tensor_I a = A.double_scalar_product(LeviCivita()) * 0.5;
    cout << "a = " << endl;
    a.show(true);
    return 0;
}