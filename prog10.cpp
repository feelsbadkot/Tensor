/*
ЗАДАЧА 10:
Дан базис 
e1 = i;
e2 = 2j;
e3 = j + k;
и контра-ковариантные компоненты тензора второго ранга
1 7 3
2 9 0
3 -4 4
Найти обратный тензор, если это возможно.
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
    Tensor_II T = Tensor_II(e, 1, comps);
    Tensor_II I = Tensor_II(e, 0, T.e().inv_g_matrix());
    // Решение
    //Tensor_II T_ = Tensor_II(e, 1, inverse_matrix(3, T.contra_co_comps()));
    //T_.show(1);
    //cout << endl;
    //(T_ * T).show(0);
    //cout << endl;
    //(T * T_).show(0);
    //cout << endl;
    double J1 = T.trace();
    double J2 = 0.5 * (pow(T.trace(), 2) - T.pow(2).trace());
    double J3 = (1 / 3.0) * T.pow(3).trace() - 0.5 * T.trace() * T.pow(2).trace() + (1 / 6.0) * pow(T.trace(), 3);
    cout << "J1 = " << J1 << endl;
    cout << "J2 = " << J2 << endl;
    cout << "J3 = " << J3 << endl;
    Tensor_II T_ = (T.pow(2) - T * J1 + I * J2) * (1 / J3);
    cout << "T^(-1) = " << endl;
    T_.show(1);
    I.show(3);
    (T * T_).show(3);
    return 0;  
}