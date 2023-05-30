/*
ЗАДАЧА 6:
Дан базис 
e1 = i;
e2 = 2j;
e3 = j + k;
Даны дважды контравариантные компоненты тензора 2 ранга:
1 2 3
2 9 0
3 0 4
Найти оставшиеся типы компонент
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
    comps[0] = new double[3] {1, 2, 3};
    comps[1] = new double[3] {-4, 9, 0};
    comps[2] = new double[3] {3, 6, 4};
    Tensor_II T = Tensor_II(e, 0, comps);
    // Решение
    cout << "T^{ij}= " << endl;
    T.show(0);
    cout << endl;
    cout << "T^i_j = " << endl;
    T.show(1);
    cout << endl;
    cout << "T_i^j = " << endl;
    T.show(2);
    cout << endl;
    cout << "T_ij = " << endl;
    T.show(3);
    return 0;
}