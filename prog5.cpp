/* 
ЗАДАЧА 5:
Дан базис 
e1 = i;
e2 = 2j;
e3 = j + k;
Даны контравариантные компоненты векторов в этом базисе:
a = 6e1 - e2 - 1/2 e3;
b = e1 + e2;
c = 4e1 - 4e3;
Найти 
(a, b)
[a, b]
(a, b, c)
cos(b, c)
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
    double *a_comps = new double[3] {6, -1, -0.5};
    Tensor_I a = Tensor_I(e, true, a_comps);
    double *b_comps = new double[3] {1, 1, 0};
    Tensor_I b = Tensor_I(e, true, b_comps);
    double *c_comps = new double[3] {4, 0, -4};
    Tensor_I c = Tensor_I(e, true, c_comps);
    // Решение
    cout << "(a, b) = " << a * b << endl;
    cout << "[a, b] = ";
    a.vector_product(b).show(true);
    cout << "(a, b, c) = " << a.mixed_product(b, c) << endl;
    cout << "cos(b, c) = " << (b * c) / (b.len() * c.len());
    return 0; 
}