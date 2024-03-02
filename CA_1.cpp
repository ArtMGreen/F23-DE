#include <iostream>
#include <cmath>
#include <vector>
#include <string>


using namespace std;


double X_0 = 1;
double Y_0 = 1;
double X_f = 3;

double y(double x) {
    return exp(x) * log(x) + exp(x-1);
}

double dy_dx(double x, double y) {
    return y + exp(x)/x;
}

double max_error(vector<double> x) {
    double max_value = -1;
    for (int i = 0; i < x.size(); i++) {
        if (x[i] > max_value) max_value = x[i];
    }
    return max_value;
}

void print_vector(vector<double> vec, string vecname) {
    printf("%s=\n", vecname.c_str());
    for (int i = 0; i < vec.size() - 1; i++) {
        printf("%.5f ", vec[i]);
    }
    printf("%.5f\n", vec[vec.size()-1]);
}

void print_vector(vector<int> vec, string vecname) {
    printf("%s=\n", vecname.c_str());
    for (int i = 0; i < vec.size() - 1; i++) {
        printf("%d ", vec[i]);
    }
    printf("%d\n", vec[vec.size()-1]);
}

vector<double> generate_xi(double x_0, double x_f, int n) {
    vector<double> result;
    double h = (x_f - x_0) / (double)n;
    double x_i = x_0;
    result.push_back(x_i);
    for (int i = 1; i <= n; i++) {
        x_i += h;
        result.push_back(x_i);
    }
    return result;
}

vector<int> generate_ni(int n1, int n2) {
    vector<int> result;
    for (int i = n1; i <= n2; i++) {
        result.push_back(i);
    }
    return result;
}

vector<double> analytical_solutions(double x_0, double x_f, int n) {
    vector<double> result;
    double h = (x_f - x_0) / (double)n;
    double x_i = x_0;
    result.push_back(y(x_i));
    for (int i = 1; i <= n; i++) {
        x_i += h;
        result.push_back(y(x_i));
    }
    return result;
}

vector<double> local_errors(vector<double> model, vector<double> result) {
    vector<double> errors;
    for (int i = 0; i < model.size(); i++) {
        errors.push_back(abs(result[i] - model[i]));
    }
    return errors;
}

vector<double> euler(double x_0, double y_0, double x_f, int n) {
    vector<double> result;
    double h = (x_f - x_0) / (double)n;
    double y_i = y_0;
    double x_i = x_0;
    result.push_back(y_i);
    for (int i = 1; i <= n; i++) {
        y_i = y_i + h * dy_dx(x_i, y_i);
        x_i += h;
        result.push_back(y_i);
    }
    return result;
}

vector<double> improved_euler(double x_0, double y_0, double x_f, int n) {
    vector<double> result;
    double h = (x_f - x_0) / (double)n;
    double y_i = y_0;
    double x_i = x_0;
    result.push_back(y_i);
    for (int i = 1; i <= n; i++) {
        double k1i = dy_dx(x_i, y_i);
        double k2i = dy_dx(x_i + h, y_i + h * k1i);
        y_i = y_i + (h/2) * (k1i + k2i);
        x_i += h;
        result.push_back(y_i);
    }
    return result;
}

vector<double> runge_kutta_4th_order(double x_0, double y_0, double x_f, int n) {
    vector<double> result;
    double h = (x_f - x_0) / (double)n;
    double y_i = y_0;
    double x_i = x_0;
    result.push_back(y_i);
    for (int i = 1; i <= n; i++) {
        double k1i = dy_dx(x_i, y_i);
        double k2i = dy_dx(x_i + (h/2), y_i + (h/2) * k1i);
        double k3i = dy_dx(x_i + (h/2), y_i + (h/2) * k2i);
        double k4i = dy_dx(x_i + h, y_i + h * k3i);
        y_i = y_i + (h/6) * (k1i + 2 * k2i + 2 * k3i + k4i);
        x_i += h;
        result.push_back(y_i);
    }
    return result;
}

int CA_1() {
    int n, n1, n2;
    cin >> n;
    cin >> n1 >> n2;
    int task;
    cin >> task;
    n -= 1;


    vector<double> xi = generate_xi(X_0, X_f, n);
    if (task == 1) {
        vector<double> y_xi = analytical_solutions(X_0, X_f, n);
        print_vector(xi, "xi");
        print_vector(y_xi, "y(xi)");
    } else if (task == 2) {
        vector<double> euler_yi = euler(X_0, Y_0, X_f, n);
        print_vector(xi, "xi");
        print_vector(euler_yi, "Euler_yi");
    }
    else if (task == 3) {
        vector<double> improved_euler_yi = improved_euler(X_0, Y_0, X_f, n);
        print_vector(xi, "xi");
        print_vector(improved_euler_yi, "iEuler_yi");
    }
    else if (task == 4) {
        vector<double> rk4 = runge_kutta_4th_order(X_0, Y_0, X_f, n);
        print_vector(xi, "xi");
        print_vector(rk4, "RK4_yi");
    }
    else if (task == 5) {
        vector<double> y_xi = analytical_solutions(X_0, X_f, n);
        vector<double> euler_yi = euler(X_0, Y_0, X_f, n);
        vector<double> euler_LE_xi = local_errors(y_xi, euler_yi);
        print_vector(xi, "xi");
        print_vector(euler_LE_xi, "Euler_LE(xi)");
    }
    else if (task == 6) {
        vector<double> y_xi = analytical_solutions(X_0, X_f, n);
        vector<double> improved_euler_yi = improved_euler(X_0, Y_0, X_f, n);
        vector<double> ieuler_LE_xi = local_errors(y_xi, improved_euler_yi);
        print_vector(xi, "xi");
        print_vector(ieuler_LE_xi, "iEuler_LE(xi)");
    }
    else if (task == 7) {
        vector<double> y_xi = analytical_solutions(X_0, X_f, n);
        vector<double> rk4 = runge_kutta_4th_order(X_0, Y_0, X_f, n);
        vector<double> rk4_LE_xi = local_errors(y_xi, rk4);
        print_vector(xi, "xi");
        print_vector(rk4_LE_xi, "RK4_LE(xi)");
    }
    else if (task == 8) {
        vector<double> GE;
        for (int n_i = n1; n_i <= n2; n_i++) {
            vector<double> y_xi = analytical_solutions(X_0, X_f, n_i-1);
            vector<double> euler_yi = euler(X_0, Y_0, X_f, n_i-1);
            vector<double> euler_LE_xi = local_errors(y_xi, euler_yi);
            GE.push_back(max_error(euler_LE_xi));
        }
        vector<int> ni = generate_ni(n1, n2);
        print_vector(ni, "ni");
        print_vector(GE, "Euler_GE(n)");
    }
    else if (task == 9) {
        vector<double> GE;
        for (int n_i = n1; n_i <= n2; n_i++) {
            vector<double> y_xi = analytical_solutions(X_0, X_f, n_i-1);
            vector<double> improved_euler_yi = improved_euler(X_0, Y_0, X_f, n_i-1);
            vector<double> ieuler_LE_xi = local_errors(y_xi, improved_euler_yi);
            GE.push_back(max_error(ieuler_LE_xi));
        }
        vector<int> ni = generate_ni(n1, n2);
        print_vector(ni, "ni");
        print_vector(GE, "iEuler_GE(n)");
    }
    else if (task == 10) {
        vector<double> GE;
        for (int n_i = n1; n_i <= n2; n_i++) {
            vector<double> y_xi = analytical_solutions(X_0, X_f, n_i-1);
            vector<double> rk4_yi = runge_kutta_4th_order(X_0, Y_0, X_f, n_i-1);
            vector<double> rk4_LE_xi = local_errors(y_xi, rk4_yi);
            GE.push_back(max_error(rk4_LE_xi));
        }
        vector<int> ni = generate_ni(n1, n2);
        print_vector(ni, "ni");
        print_vector(GE, "RK4_GE(n)");
    }
    return 0;
}
