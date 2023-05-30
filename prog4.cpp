/*
ЗАДАЧА 4:
Дан базис 
e1 = i;
e2 = 2j;
e3 = j + k;
Даны контравариантные компоненты вектора в этом базисе: 
x = 5e1 - 1/2 e2 - 3e3;
Дан новый базис 
e1' = e1;
e2' = e1 + e2;
e3' = e1 + e2 + e3;
Найти компоненты вектора x в новом базисе
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
    double *comps = new double[3] {5, -0.5, -3};
    Tensor_I x = Tensor_I(e, true, comps);
    cout << "x = ";
    x.show(true);
    x.show(false);
    cout << "|x| = " << x.len() << endl;
    Basis e_new = Basis(e1, e1 + e2, e1 + e2 + e3);
    cout << "e'1 = ";
    e_new.e1().show();
    cout << "e'2 = ";
    e_new.e2().show();
    cout << "e'3 = ";
    e_new.e3().show();
    // Строим матрицу перехода:
    double **a = new double*[3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i] = new double[3]{0};
        }
    }
    a[0][0] = 1;
    a[0][1] = 1;
    a[0][2] = 1;
    a[1][1] = 1;
    a[1][2] = 1;
    a[2][2] = 1;
    double **b = inverse_matrix(3, a);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << b[i][j] << " ";
        }
        cout << endl;
    }
    // Преобразование компонент по известному закону: 
    double *comps_new = new double[3]{0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            comps_new[i] += b[i][j] * x.contra_comps()[j];
        }
    }
    Tensor_I x_new = Tensor_I(e_new, true, comps_new);
    cout << "x_new = ";
    x_new.show(true);
    x_new.show(false);
    cout  << "|x| = " << x_new.len() << endl;
}