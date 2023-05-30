/*
ЗАДАЧА 2:
Дан базис: 
e1 = i;
e2 = 2j;
e3 = j + k;
Даны контравариантные компоненты вектора:
a = 3e1 + 5e2 - e3;
Найти ковариантные компоненты вектора и его модуль
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
    double *comps = new double[3] {3, 5, -1};
    Tensor_I a = Tensor_I(e, true, comps);
    cout << "a = ";
    a.show(true);
    // Ищем ковариантные компоненты
    cout << "a = ";
    a.show(false);
    // Ищем длину
    cout << "|a| = " << a.len() << endl; 
    return 0;
}