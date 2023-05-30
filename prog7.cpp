/*
ЗАДАЧА 7:
Дан базис 
e1 = i;
e2 = 2j;
e3 = j + k;
Даны дважды контравариантные компоненты тензора 2 ранга:
1 7 3
2 9 0
3 -4 4
Найти его симметричную и антисимметричную части
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
    double **comps = new double*[3];
    comps[0] = new double[3] {1, 7, 3};
    comps[1] = new double[3] {2, 9, 0};
    comps[2] = new double[3] {3, -4, 4};
    Tensor_II T = Tensor_II(e, 3, comps);
    // Решение
    Tensor_II Ts = T.symmetric_part();
    cout << "Ts = " << endl;
    Ts.show(0);
    Tensor_II Ta = T.antisymmetric_part();
    cout << "Ta = " << endl;
    Ta.show(0);
    return 0;
}