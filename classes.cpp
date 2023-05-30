#include <iostream>
#include <cmath>

// функция возвразающая матрицу без i-ой строки и j-ого столбца 
double **get_minor_matrix(int n, double **A, int i, int j) {
    // выделяем память под новую уменьшенную матрицу
    double **minor_matr = new double*[n - 1];
    for (int k = 0; k < n - 1; k++) {
        minor_matr[k] = new double[n - 1];
    }
    // счетчики для новой матрицы
    int k = 0, m = 0;
    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++) {
            // записываем в новую матрицу элементы старой за исключением i-ой строки и j-ого столбца
            if (r != i && c != j) {
                minor_matr[k][m] = A[r][c];
                m++;
                // когда место в строке закончилось, начинаем заполнять новую
                if (m == n - 1) {
                    m = 0;
                    k++;
                } 
            }
        }
    }
    return minor_matr;
}

// функция, рекурсивно вычисляющая определитель
double determinant(int n, double **A) {
    // база рекурсии
    double det = 0;
    if (n == 1) {
        return A[0][0];
    }
    // раскладываем матрицу по первому столбцу
    for (int i = 0; i < n; i++) {
        double **minor_matr = get_minor_matrix(n, A, i, 1);
        det += (pow(-1, i + 1) * A[i][1] * determinant(n - 1, minor_matr));
        for (int k = 0; k < n - 1; k++)
            delete minor_matr[k];
        delete minor_matr;
    }
    return det;
}

double **inverse_matrix(int n, double **A) {
    // выделяем память под матрицу алгебраических дополнений 
    double **alg_add_A = new double*[n];
    for (int i = 0; i < n; i++) {
        alg_add_A[i] = new double[n];
    }
    // вычисляем определитель заданной матрицы
    double det_A = determinant(n, A);
    // составляем матрицу дополнений
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double **minor_matr = get_minor_matrix(n, A, i, j);
            alg_add_A[i][j] = pow(-1, i + j) * determinant(n - 1, minor_matr);
            for (int k = 0; k < n - 1; k++) {
                delete minor_matr[k];
            }
            delete minor_matr;
        }
    }
    // выделяем память под обратную матрицу
    double **inv_A = new double*[n];
    for (int i = 0; i < n; i++) {
        inv_A[i] = new double[n];
    }
    // через транспонирование и поэлементное деление на определитель вычисляем обратную матрицу
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv_A[i][j] = alg_add_A[j][i] / det_A;
        }
    }
    for (int i = 0; i < n; i++) {
        delete alg_add_A[i];
    }
    delete alg_add_A;
    return inv_A;
}

class Vector {
    private:
        double x_;
        double y_;
        double z_;
    public:
        Vector() = default;
        Vector(double, double, double);
        ~Vector() = default;
        double x();
        double y();
        double z();
        Vector operator+(Vector);
        Vector operator-(Vector);
        void operator+=(Vector);
        void operator-=(Vector);
        bool operator==(Vector);
        Vector operator*(double);
        double operator*(Vector);
        Vector operator/(Vector);
        double mixed_product(Vector, Vector);
        double len();
        void show();
        bool is_zero();
};

class Basis {
    private: 
        Vector e1_;
        Vector e2_;
        Vector e3_;
        Vector *e_;
        double **g_matrix_;
        double **inv_g_matrix_;
        double g_;
        double inv_g_;

    public:
        Basis() = default; 
        Basis(Vector, Vector, Vector);
        double **g_matrix();
        double **inv_g_matrix();
        double g();
        double inv_g();
        Vector e1();
        Vector e2();
        Vector e3();
        Vector *e();
        bool operator==(Basis);
        ~Basis() = default;
};

class Tensor_I {
    private:
        const short rank = 1;
        Basis basis_;
        double *co_comps_;
        double *contra_comps_;
    
    public:
        Tensor_I() = default;
        Tensor_I(Basis, bool, double*);
        ~Tensor_I() = default;
        double *co_comps();
        double *contra_comps();
        Tensor_I operator+(Tensor_I);
        Tensor_I operator-(Tensor_I);
        void operator+=(Tensor_I);
        void operator-=(Tensor_I);
        bool operator==(Tensor_I);
        // умножение на число
        template<class Type>
        Tensor_I operator*(Type);
        // скалярное умножение
        double operator*(Tensor_I);
        Tensor_I vector_product(Tensor_I);
        double mixed_product(Tensor_I, Tensor_I);
        double len();
        void show(bool);
        bool is_zero();
};

class LeviCivita {
    private:
        int ***comps_;
    public:
        LeviCivita();
        int ***comps();
        int comps(int, int, int);
};

class Tensor_II {
    private:
        const short rank = 2;
        Basis basis_;
        double **double_contra_comps_;
        double **double_co_comps_;
        double **co_contra_comps_;
        double **contra_co_comps_;
    public:
        Tensor_II() = default;
        ~Tensor_II() = default;
        Tensor_II(Basis, int, double**);
        double **double_contra_comps();
        double **double_co_comps();
        double **co_contra_comps();
        double **contra_co_comps();
        Basis e();
        Tensor_II operator+(Tensor_II);
        Tensor_II operator-(Tensor_II);
        void operator+=(Tensor_II);
        void operator-=(Tensor_II);
        bool operator==(Tensor_II); 
        template<class Type>
        Tensor_II operator*(Type);
        Tensor_I operator*(Tensor_I);
        Tensor_II operator*(Tensor_II);
        Tensor_II transpose();
        Tensor_II pow(int n);
        Tensor_II symmetric_part();
        Tensor_II antisymmetric_part();
        Tensor_II sphrerical_part();
        Tensor_II deviatoric_part();
        double trace();
        double total_scalar_product(Tensor_II);
        double double_scalar_product(Tensor_II);
        Tensor_I double_scalar_product(LeviCivita);
        void show(int);
        bool is_symmetric();
        bool is_antisymmetric();
};

Vector::Vector(double x, double y, double z): x_(x), y_(y), z_(z) {}

double Vector::x() {
    return this->x_;
}

double Vector::y() {
    return this->y_;
}

double Vector::z() {
    return this->z_;
}

Vector Vector::operator+(Vector other) {
    return Vector(this->x_ + other.x_, this->y_ + other.y_, this->z_ + other.z_);
}

Vector Vector::operator-(Vector other) {
    return Vector(this->x_ - other.x_, this->y_ - other.y_, this->z_ - other.z_);
}

void Vector::operator+=(Vector other) {
    this->x_ += other.x_;
    this->y_ += other.y_;
    this->z_ += other.z_;
}

void Vector::operator-=(Vector other) {
    this->x_ -= other.x_;
    this->y_ -= other.y_;
    this->z_ -= other.z_;
}

bool Vector::operator==(Vector other) {
    return this->x_ == other.x_ && this->y_ == other.y_ && this->z_ == other.z_;
}

Vector Vector::operator*(double num) {
    return Vector(num * this->x_, num * this->y_, num * this->z_);
}

double Vector::operator*(Vector other) {
    return this->x_ * other.x_ + this->y_ * other.y_ + this->z_ * other.z_;
}

double Vector::len() {
    return sqrt((*this) * (*this));
}

Vector Vector::operator/(Vector other) {
    double x = (this->y_ * other.z_) - (other.y_ * this->z_);
    double y = -((this->x_ - other.z_) - (other.x_ - this->z_));
    double z = (this->x_ * other.y_) - (other.x_ * this->y_);
    return Vector(x, y, z);
}

double Vector::mixed_product(Vector second, Vector third) {
    return (*this) * (second / third);
}

void Vector::show() {
    std::cout << "{" << this->x() << "; " << this->y() << "; " << this->z() << "}" << std::endl; 
}

bool Vector::is_zero() {
    if (this->x_ == 0 && this->y_ == 0 && this->z_ == 0) {
        return true;
    }  
    return false;
}

Basis::Basis(Vector e1, Vector e2, Vector e3) {
    if (e1.mixed_product(e2, e3)) {
        e1_ = e1;
        e2_ = e2;
        e3_ = e3;
        e_ = new Vector[3] {e1_, e2_, e3_};
        g_matrix_ = new double*[3];
        for (int i = 0; i < 3; i++) {
            g_matrix_[i] = new double[3];
            for (int j = 0; j < 3; j++) {
                g_matrix_[i][j] = e_[i] * e_[j];
            }
        }
        g_ = determinant(3, g_matrix_);
        inv_g_matrix_ = inverse_matrix(3, g_matrix_);
        inv_g_ = 1 / g_;
    }
    else {
        std::cout << "It's impossible to build a basis on {e1, e2, e3}" << std::endl;
        delete this;
    }
}

double** Basis::g_matrix() {
    return this->g_matrix_;
}

double** Basis::inv_g_matrix() {
    return this->inv_g_matrix_;
}

double Basis::g() {
    return this->g_;
}

double Basis::inv_g() {
    return this->inv_g_;
}

Vector Basis::e1() {
    return this->e1_;
}

Vector Basis::e2() {
    return this->e2_;
}

Vector Basis::e3() {
    return this->e3_;
}

Vector* Basis::e() {
    return this->e_;
}

bool Basis::operator==(Basis other) {
    return (this->e1_ == other.e1_ && this->e2_ == other.e2_ && this->e3_ == other.e3_);
}

Tensor_I::Tensor_I(Basis basis, bool is_main_basis, double *comps) {
    this->co_comps_ = new double[3]{0};
    this->contra_comps_ = new double[3]{0};
    this->basis_ = basis;
    if (is_main_basis) { // значит даны контравариантные компоненты
        for (int i = 0; i < 3; i++) {
            this->contra_comps_[i] = comps[i];
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                this->co_comps_[i] += this->basis_.g_matrix()[i][j] * this->contra_comps_[j];
            }
        }
    }
    else { // значит даны ковариантные компоненты
        for (int i = 0; i < 3; i++) {
            this->co_comps_[i] = comps[i];
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                this->contra_comps_[i] += this->basis_.inv_g_matrix()[i][j] * this->co_comps_[j];
            }
        }
    }
}

double* Tensor_I::co_comps() {
    return this->co_comps_;
}

double* Tensor_I::contra_comps() {
    return this->contra_comps_;
}

Tensor_I Tensor_I::operator+(Tensor_I other) {
    if (this->basis_ == other.basis_) {
        double *contra_comps = new double[3];
        for (int i = 0; i < 3; i++) {
            contra_comps[i] = this->contra_comps_[i] + other.contra_comps_[i];
        }
        return Tensor_I(this->basis_, true, contra_comps);
    }
    return Tensor_I();
}

Tensor_I Tensor_I::operator-(Tensor_I other) {
    if (this->basis_ == other.basis_) {
        double *contra_comps = new double[3];
        for (int i = 0; i < 3; i++) {
            contra_comps[i] = this->contra_comps_[i] - other.contra_comps_[i];
        }
        return Tensor_I(this->basis_, true, contra_comps);
    }
    return Tensor_I();
}

void Tensor_I::operator+=(Tensor_I other) {
    if (this->basis_ == other.basis_) {
        for (int i = 0; i < 3; i++) {
            this->contra_comps_[i] += other.contra_comps_[i];
            this->co_comps_[i] += other.co_comps_[i];
        }
    } 
}

void Tensor_I::operator-=(Tensor_I other) {
    if (this->basis_ == other.basis_) {
        for (int i = 0; i < 3; i++) {
            this->contra_comps_[i] -= other.contra_comps_[i];
            this->co_comps_[i] -= other.co_comps_[i];
        }
    } 
}

bool Tensor_I::operator==(Tensor_I other) {
    if (this->basis_ == other.basis_) {
        for (int i = 0; i < 3; i++) {
            if (this->co_comps_[i] != other.co_comps_[i] || this->contra_comps_[i] != other.contra_comps_[i]) {
                return false;
            }
        }
        return true;
    }
    return false;
}

// умножение на число
template<class Type>
Tensor_I Tensor_I::operator*(Type num) {
    double *comps = new double[3] {0};
    comps[0] = num * this->contra_comps_[0];
    comps[1] = num * this->contra_comps_[1];
    comps[2] = num * this->contra_comps_[2];
    //std::cout << " " << comps[0] << " " << comps[1] << " " << comps[2] << std::endl;
    return Tensor_I(this->basis_, true, comps);
}

// скалярное умножение
double Tensor_I::operator*(Tensor_I other) {
    double result = 0;
    for (int i = 0; i < 3; i++) {
        result += (this->contra_comps_[i] * other.co_comps_[i]);
    }
    return result;
}

//Tensor_I Tensor_I::operator*(Tensor_II other) {
//    double *comps = new double[3]{0};
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < 3; j++) {
//            for (int k = 0; k < 3; k++) {
//                comps[k] += this->contra_comps_[i] * other.double_contra_comps()[j][k] * this->basis_.g_matrix()[i][j];
//            }
//        }
//    }
//    return Tensor_I(this->basis_, true, comps);
//}

double Tensor_I::len() {
    return sqrt((*this) * (*this));
}

void Tensor_I::show(bool is_main_basis) {
    if (is_main_basis) {
        std::cout << "{" << contra_comps_[0] << "; " << contra_comps_[1] << "; " << contra_comps_[2] << "}" << std::endl;
    }
    else {
        std::cout << "{" << co_comps_[0] << "; " << co_comps_[1] << "; " << co_comps_[2] << "}" << std::endl;
    }
}

Tensor_I Tensor_I::vector_product(Tensor_I other) {     
    if (this->basis_ == other.basis_) {
        double v1 = sqrt(this->basis_.g()) * ((this->contra_comps_[1] * other.contra_comps_[2]) - (other.contra_comps_[1] * this->contra_comps_[2]));
        double v2 = -sqrt(this->basis_.g()) * ((this->contra_comps_[0] * other.contra_comps_[2]) - (other.contra_comps_[0] * this->contra_comps_[2]));
        double v3 = sqrt(this->basis_.g()) * ((this->contra_comps_[0] * other.contra_comps_[1]) - (other.contra_comps_[0] * this->contra_comps_[1]));
        double *comps = new double[3] {v1, v2, v3};
        return Tensor_I(this->basis_, true, comps);
    }
    std::cout << "Basis error" << std::endl;
    return Tensor_I();
}

double Tensor_I::mixed_product(Tensor_I second, Tensor_I third) {
    if (this->basis_ == second.basis_ &&  second.basis_ == third.basis_) {
        return (*this) * second.vector_product(third);
    }
    std::cout << "Basis error" << std::endl;
    return INT_MAX;
}

bool Tensor_I::is_zero() {
    if (this->contra_comps_[0] == 0 && this->contra_comps_[1] == 0 && this->contra_comps_[2] == 0) {
        return true;
    }
    return false;
}

LeviCivita::LeviCivita() {
    comps_ = new int**[3];
    for (int i = 0; i < 3; i++) {
        comps_[i] = new int*[3];
        for (int j = 0; j < 3; j++) {
            comps_[i][j] = new int[3]{0};
            for (int k = 0; k < 3; k++) {
                if ((i == 0 && j == 1 && k == 2) || (i == 2 && j == 0 && k == 1) || (i == 1 && j == 2 && k == 0)) {
                    comps_[i][j][k] = 1;
                } 
                if ((i == 0 && j == 2 && k == 1) || (i == 2 && j == 1 && k == 0) || (i == 1 && j == 0 && k == 2)) {
                    comps_[i][j][k] = -1;
                }
            }
        }
    }
}

int*** LeviCivita::comps() {
    return this->comps_;
}

int LeviCivita::comps(int i, int j, int k) {
    return this->comps_[i][j][k];
}

Tensor_II::Tensor_II(Basis basis, int type_of_comps, double **comps) {
    this->basis_ = basis;
    // Считаем компоненты
    this->double_contra_comps_ = new double*[3];
    this->double_co_comps_ = new double*[3];
    this->co_contra_comps_ = new double*[3];
    this->contra_co_comps_ = new double*[3];
    double **delta_cr = new double*[3];
    for (int i = 0; i < 3; i++) {
        this->double_contra_comps_[i] = new double[3]{0};
        this->double_co_comps_[i] = new double[3]{0};
        this->co_contra_comps_[i] = new double[3]{0};
        this->contra_co_comps_[i] = new double[3]{0};
        delta_cr[i] = new double[3]{0};
        delta_cr[i][i] = 1;
    }
    // дважды контравариантные 
    if (type_of_comps == 0) {
        double_contra_comps_ = comps;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    for (int m = 0; m < 3; m++) {
                        this->contra_co_comps_[k][m] += this->basis_.g_matrix()[j][m] * delta_cr[i][k] * this->double_contra_comps_[i][j];
                        this->co_contra_comps_[k][m] += this->basis_.g_matrix()[i][k] * delta_cr[j][m] * this->double_contra_comps_[i][j];
                        this->double_co_comps_[k][m] += this->basis_.g_matrix()[i][k] * this->basis_.g_matrix()[j][m] * this->double_contra_comps_[i][j];
                    }
                }
            }
        }
    }
    // контра-ковариантные  
    if (type_of_comps == 1) {
        contra_co_comps_ = comps;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    for (int m = 0; m < 3; m++) {
                        double_co_comps_[k][m] += basis_.g_matrix()[i][k] * delta_cr[j][m] * contra_co_comps_[i][j]; 
                        double_contra_comps_[k][m] += basis_.inv_g_matrix()[j][m] * delta_cr[i][k] * contra_co_comps_[i][j];
                        co_contra_comps_[k][m] += basis_.g_matrix()[i][k] * basis_.inv_g_matrix()[j][m] * contra_co_comps_[i][j];                
                    }
                }
            }
        }
    }
    // ко-контравариантные
    if (type_of_comps == 2) {
        co_contra_comps_ = comps;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    for (int m = 0; m < 3; m++) {
                        double_contra_comps_[k][m] += basis_.inv_g_matrix()[i][k] * delta_cr[j][m] * co_contra_comps_[i][j];
                        double_co_comps_[k][m] += basis_.g_matrix()[j][m] * delta_cr[i][k] * co_contra_comps_[i][j];
                        contra_co_comps_[k][m] += basis_.inv_g_matrix()[i][k] * basis_.g_matrix()[j][m] * co_contra_comps_[i][j];
                    }
                }
            }
        }
    }
    // дважды ковариантные компоненты
    if (type_of_comps == 3) {
        double_co_comps_ = comps;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    for (int m = 0; m < 3; m++) {
                        co_contra_comps_[k][m] += basis_.inv_g_matrix()[i][k] * delta_cr[j][m] * double_co_comps_[i][j];
                        contra_co_comps_[k][m] += basis_.inv_g_matrix()[j][m] * delta_cr[i][k] * double_co_comps_[i][j];
                        double_contra_comps_[k][m] += basis_.inv_g_matrix()[i][k] * basis_.inv_g_matrix()[j][m] * double_co_comps_[i][j]; 
                    }
                }
            }
        }
    }
}

double** Tensor_II::double_contra_comps() {
    return this->double_contra_comps_;
}

double** Tensor_II::double_co_comps() {
    return this->double_co_comps_;
}

double** Tensor_II::co_contra_comps() {
    return this->co_contra_comps_;
}

double** Tensor_II::contra_co_comps() {
    return this->contra_co_comps_;
}

Basis Tensor_II::e() {
    return this->basis_;
}

Tensor_II Tensor_II::operator+(Tensor_II other) {
    if (this->basis_ == other.basis_) {
        double **comps = new double*[3];
        for (int i = 0; i < 3; i++) {
            comps[i] = new double[3]{0};
            for (int j = 0; j < 3; j++) {
                comps[i][j] = this->double_contra_comps_[i][j] + other.double_contra_comps_[i][j];
            }
        }
        return Tensor_II(this->basis_, 0, comps);
    } 
    return Tensor_II();
}

Tensor_II Tensor_II::operator-(Tensor_II other) {
    if (this->basis_ == other.basis_) {
        double **comps = new double*[3];
        for (int i = 0; i < 3; i++) {
            comps[i] = new double[3]{0};
            for (int j = 0; j < 3; j++) {
                comps[i][j] = this->double_contra_comps_[i][j] - other.double_contra_comps_[i][j];
            }
        }
        return Tensor_II(this->basis_, 0, comps);
    }
    return Tensor_II();
}

void Tensor_II::operator+=(Tensor_II other) {
    if (this->basis_ == other.basis_) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                this->double_contra_comps_[i][j] += other.double_contra_comps_[i][j];
                this->double_co_comps_[i][j] += other.double_co_comps_[i][j];
                this->co_contra_comps_[i][j] += other.co_contra_comps_[i][j];
                this->contra_co_comps_[i][j] += other.contra_co_comps_[i][j];
            }
        }
    }
}

void Tensor_II::operator-=(Tensor_II other) {
    if (this->basis_ == other.basis_) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                this->double_contra_comps_[i][j] -= other.double_contra_comps_[i][j];
                this->double_co_comps_[i][j] -= other.double_co_comps_[i][j];
                this->co_contra_comps_[i][j] -= other.co_contra_comps_[i][j];
                this->contra_co_comps_[i][j] -= other.contra_co_comps_[i][j];
            }
        }
    }
}

bool Tensor_II::operator==(Tensor_II other) {
    if (this->basis_ == other.basis_) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (this->double_contra_comps_[i][j] != other.double_contra_comps_[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }
    return false;
}

template<class Type>
Tensor_II Tensor_II::operator*(Type num) {
    double **comps = new double*[3];
    for (int i = 0; i < 3; i++) {
        comps[i] = new double[3]{0};
        for (int j = 0; j < 3; j++) {
            comps[i][j] = num * this->double_contra_comps_[i][j];
        }
    }
    return Tensor_II(this->basis_, 0, comps);
}

Tensor_I Tensor_II::operator*(Tensor_I other) {
    double *comps = new double[3]{0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                comps[i] += this->double_contra_comps_[i][j] * other.contra_comps()[k] * this->basis_.g_matrix()[j][k];
            }
        }
    }
    return Tensor_I(this->basis_, true, comps);
}

Tensor_II Tensor_II::operator*(Tensor_II other) {
    if (this->basis_ == other.basis_) {
        double **comps = new double*[3];
        for (int i = 0; i < 3; i++) {
            comps[i] = new double[3]{0};
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    for (int m = 0; m < 3; m++) {
                        comps[i][m] += this->double_contra_comps_[i][j] * other.double_contra_comps_[k][m] * this->basis_.g_matrix()[j][k];
                    }
                }
            }
        }
        return Tensor_II(this->basis_, 0, comps);
    }
    return Tensor_II();
} 

Tensor_II Tensor_II::transpose() {
    double **comps = new double*[3];
    for (int i = 0; i < 3; i++) {
        comps[i] = new double[3]{0};
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            comps[i][j] = this->double_contra_comps_[j][i];
        }
    }
    return Tensor_II(this->basis_, 0, comps);
}

Tensor_II Tensor_II::pow(int n) {
    if (n > 0) {
        double **comps = new double*[3];
        for (int i = 0; i < 3; i++) {
            comps[i] = new double[3]{0};
        }
        if (n == 1) {
            return (*this);
        }
        if (n == 2) {
            (*this) * (*this);
        }
        return (*this) * this->pow(n - 1);
    }
    return Tensor_II();
}

Tensor_II Tensor_II::symmetric_part() {
    return ((*this) + this->transpose()) * 0.5;
}

Tensor_II Tensor_II::antisymmetric_part() {
    return ((*this) - this->transpose()) * 0.5;
}

double Tensor_II::trace() {
    double trace = 0;
    for (int i = 0; i < 3; i++) {
        trace += contra_co_comps_[i][i];
    }
    return trace;
}

Tensor_II Tensor_II::sphrerical_part() {
    double **comps = new double*[3];
    for (int i = 0; i < 3; i++) {
        comps[i] = new double[3]{0};
        comps[i][i] = 1;
    }
    return Tensor_II(this->basis_, 1, comps) * this->trace() * (1 / 3.0);
}

Tensor_II Tensor_II::deviatoric_part() {   
    return (*this) - this->sphrerical_part();
}

double Tensor_II::total_scalar_product(Tensor_II other) {
    double result = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int m = 0; m < 3; m++) {
                    result += this->double_contra_comps_[i][j] * other.double_contra_comps_[k][m] * this->basis_.g_matrix()[i][k] * this->basis_.g_matrix()[j][m];
                }
            }
        }
    }
    return result;
}

double Tensor_II::double_scalar_product(Tensor_II other) {
    double result = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int m = 0; m < 3; m++) {
                    result += this->double_contra_comps_[i][j] * other.double_contra_comps_[k][m] * this->basis_.g_matrix()[i][m] * this->basis_.g_matrix()[j][k];
                }
            }
        }
    }
    return result;
}

bool Tensor_II::is_symmetric() {
    if ((*this) == this->transpose()) {
        return true;
    }
    return false;
}

bool Tensor_II::is_antisymmetric() {
    if ((*this) == this->transpose() * (-1)) {
        return true;
    }
    return false;
}

Tensor_I Tensor_II::double_scalar_product(LeviCivita E) {
    double *comps = new double[3]{0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int m = 0; m < 3; m++) {
                    for (int l = 0; l < 3; l++) {
                        comps[l] += this->double_contra_comps_[i][j] * E.comps(k, m, l) * this->basis_.g_matrix()[i][m] * this->basis_.g_matrix()[j][k]; 
                    }
                }
            } 
        }
    }
    return Tensor_I(this->basis_, true, comps);
}   

void Tensor_II::show(int type_of_comps) {
    if (type_of_comps == 0) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << this->double_contra_comps_[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    if (type_of_comps == 1) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << this->contra_co_comps_[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    if (type_of_comps == 2) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << this->co_contra_comps_[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    if (type_of_comps == 3) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << this->double_co_comps_[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
}