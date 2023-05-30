/*
ЗАДАЧА 1:
Дан базис: 
e1 = i;
e2 = 2j;
e3 = j + k;
Найти взаимный базис
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
    // Фундаментальная матрица: 
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << e.g_matrix()[i][j] << " ";
        }
        cout << endl;
    }
    // Ищем взаимный базис:
    Vector inv_e1 = Vector(0, 0, 0);
    Vector inv_e2 = Vector(0, 0, 0);
    Vector inv_e3 = Vector(0, 0, 0);
    Vector inv_e[3] = {inv_e1, inv_e2, inv_e3};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            inv_e[i] += e.e()[j] * e.inv_g_matrix()[i][j];
        }
        cout << "dual e" << i + 1 << " = ";
        inv_e[i].show();
    }
    Basis dual_e = Basis(inv_e[0], inv_e[1], inv_e[2]);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << dual_e.g_matrix()[i][j] << " ";
        }
        cout << endl;
    }
    return 0;
}