#include <iostream>
#include <cmath>
#include <vector>
#include <string>


using namespace std;


struct vector_pair {
    vector<double> a;
    vector<double> b;
};


double X_0 = -2;
double Y_0 = -2;
double T_0 = 0;
double T_f = 0.55;

double dy_dt(double x, double y, double t) {
    return 2*x*x - 2*x*y;
}

double dx_dt(double x, double y, double t) {
    return 2 + y - x*x;
}

void print_vector(vector<double> vec, string vecname) {
    printf("%s=\n", vecname.c_str());
    for (int i = 0; i < vec.size() - 1; i++) {
        printf("%.5f, ", vec[i]);
    }
    printf("%.5f\n", vec[vec.size()-1]);
}

vector<double> generate_ti(double t_0, double t_f, int n) {
    vector<double> result;
    double h = (t_f - t_0) / (double)n;
    double t_i = t_0;
    result.push_back(t_i);
    for (int i = 1; i <= n; i++) {
        t_i += h;
        result.push_back(t_i);
    }
    return result;
}

vector_pair euler(double x_0, double y_0, double t_0, double t_f, int n) {
    vector<double> result_x;
    vector<double> result_y;
    double h = (t_f - t_0) / (double)n;
    double y_i = y_0;
    double x_i = x_0;
    double t_i = t_0;
    result_x.push_back(x_i);
    result_y.push_back(y_i);
    for (int i = 1; i <= n; i++) {
        double y_next = y_i + h * dy_dt(x_i, y_i, t_i);
        double x_next = x_i + h * dx_dt(x_i, y_i, t_i);
        t_i += h;

        y_i = y_next;
        x_i = x_next;
        result_x.push_back(x_i);
        result_y.push_back(y_i);
    }
    vector_pair result = vector_pair{result_x, result_y};
    return result;
}

vector_pair improved_euler(double x_0, double y_0, double t_0, double t_f, int n) {
    vector<double> result_x;
    vector<double> result_y;
    double h = (t_f - t_0) / (double)n;
    double y_i = y_0;
    double x_i = x_0;
    double t_i = t_0;
    result_x.push_back(x_i);
    result_y.push_back(y_i);
    for (int i = 1; i <= n; i++) {
        double k1iy = dy_dt(x_i, y_i, t_i);
        double k1ix = dx_dt(x_i, y_i, t_i);

        double k2iy = dy_dt(x_i + h * k1ix, y_i + h * k1iy, t_i + h);
        double k2ix = dx_dt(x_i + h * k1ix, y_i + h * k1iy, t_i + h);

        y_i = y_i + (h/2) * (k1iy + k2iy);
        x_i = x_i + (h/2) * (k1ix + k2ix);

        t_i += h;
        result_x.push_back(x_i);
        result_y.push_back(y_i);
    }
    vector_pair result = vector_pair{result_x, result_y};
    return result;
}

vector_pair runge_kutta_4th_order(double x_0, double y_0, double t_0, double t_f, int n) {
    vector<double> result_x;
    vector<double> result_y;
    double h = (t_f - t_0) / (double)n;
    double y_i = y_0;
    double x_i = x_0;
    double t_i = t_0;
    result_x.push_back(x_i);
    result_y.push_back(y_i);
    for (int i = 1; i <= n; i++) {
        double k1iy = dy_dt(x_i, y_i, t_i);
        double k1ix = dx_dt(x_i, y_i, t_i);

        double k2iy = dy_dt(x_i + (h/2) * k1ix, y_i + (h/2) * k1iy, t_i + (h/2));
        double k2ix = dx_dt(x_i + (h/2) * k1ix, y_i + (h/2) * k1iy, t_i + (h/2));

        double k3iy = dy_dt(x_i + (h/2) * k2ix, y_i + (h/2) * k2iy, t_i + (h/2));
        double k3ix = dx_dt(x_i + (h/2) * k2ix, y_i + (h/2) * k2iy, t_i + (h/2));

        double k4iy = dy_dt(x_i + h * k3ix, y_i + h * k3iy, t_i + h);
        double k4ix = dx_dt(x_i + h * k3ix, y_i + h * k3iy, t_i + h);

        y_i = y_i + (h/6) * (k1iy + 2 * k2iy + 2 * k3iy + k4iy);
        x_i = x_i + (h/6) * (k1ix + 2 * k2ix + 2 * k3ix + k4ix);

        t_i += h;
        result_x.push_back(x_i);
        result_y.push_back(y_i);
    }
    vector_pair result = vector_pair{result_x, result_y};
    return result;
}

int CA_2() {
    int n;
    cin >> n;
    int task;
    cin >> task;
    n -= 1;

    vector<double> ti = generate_ti(T_0, T_f, n);
    if (task == 1) {
        vector_pair euler_result = euler(X_0, Y_0, T_0, T_f, n);
        vector<double> euler_xi = euler_result.a;
        vector<double> euler_yi = euler_result.b;
        print_vector(ti, "ti");
        print_vector(euler_xi, "Euler_xi");
        print_vector(euler_yi, "Euler_yi");
    } else if (task == 2) {
        vector_pair ieuler_result = improved_euler(X_0, Y_0, T_0, T_f, n);
        vector<double> ieuler_xi = ieuler_result.a;
        vector<double> ieuler_yi = ieuler_result.b;
        print_vector(ti, "ti");
        print_vector(ieuler_xi, "iEuler_xi");
        print_vector(ieuler_yi, "iEuler_yi");
    }
    else if (task == 3) {
        vector_pair rk4_result = runge_kutta_4th_order(X_0, Y_0, T_0, T_f, n);
        vector<double> rk4_xi = rk4_result.a;
        vector<double> rk4_yi = rk4_result.b;
        print_vector(ti, "ti");
        print_vector(rk4_xi, "RK4_xi");
        print_vector(rk4_yi, "RK4_yi");
    }
    return 0;
}