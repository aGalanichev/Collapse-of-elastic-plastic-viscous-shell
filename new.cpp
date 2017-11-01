#include <iostream>
#include <cmath>
#include <array>
#include <valarray>

using std::cout;
using std::endl;
using std::array;
using std::valarray;

const double PO_0 = 4530;
const double K = 123.4 * std::pow(10, 9);
const double MU = 43.4 * std::pow(10, 9);
const double a0 = std::sqrt((K+4*MU/3)/PO_0);
const double R1 = 0.02;
const double R2 = 0.03;
const double H = R2 - R1;
const double NU = 700;
const double Y = 0.71 * std::pow(10, 9);

const double P_right = -14.3 *std::pow(10, 9)/K;

const int N_x = 21;
const double h = (H)/(N_x-1);
const double t = 0.00001;

struct Point_x {
    double PO;
    double SIGMA;
    double SIGMA_R;
    double S_R;
    
    Point_x() { PO = 1; SIGMA = SIGMA_R = S_R = 0; }
};

double next_V(double v, double _r, double r, double r_, Point_x &point_l, Point_x &point_r) {
    double alpha = (point_r.PO*(r_ - r) + point_l.PO*(r - _r))/2;
    double beta  = (point_r.S_R)/((r_ + r)*point_r.PO) + (point_l.S_R)/((r + _r)*point_l.PO);

    return v + (K*t/a0/a0)*( (point_r.SIGMA_R - point_l.SIGMA_R)/alpha + 3*beta );
}

double next_V_l(double v, double r, double r_, Point_x &point, double F_l) {
    double alpha = point.PO*(r - r_);
    double beta  = (F_l - point.SIGMA)/(2*r*point.PO);

    return v + (K*t/a0/a0)*( (point.SIGMA_R - F_l)/alpha + 3*beta );
}

double next_V_r(double v, double _r, double r, Point_x &point, double F_r) {
    
    double alpha = point.PO*(r - _r);
    double beta  = (F_r - point.SIGMA)/(2*r*point.PO);
    double V = v + (K*t/a0/a0)*( (F_r - point.SIGMA)/alpha + 3*beta );
    
    return V;
}


double next_R(double r, double v) { return r + v*t; }
double get_Rc(double r_next, double r_prev) { return (r_next + r_prev)/2; }
double get_e_R(double _v, double v_, double _rc, double rc_) { return (v_ - _v)/(rc_ - _rc); }
double get_e_O(double _v, double v_, double _rc, double rc_) { return (v_ + _v)/(rc_ + _rc); }

double next_PO(Point_x &point, double e_r, double e_o) {
    return point.PO*(1 - (H*t/2/a0)*(e_r + 2*e_o))/(1 + (H*t/2/a0)*(e_r + 2*e_o));
}

double next_S_R(Point_x &point, double e_r, double e_o) {
    double S_r = point.S_R + (4*t*H*MU/(3*a0*K))*(e_r - e_o);
    double S_r_c = (S_r + point.S_R)/2;
    
    
    if(std::abs(S_r_c) - 2*Y/(3*K) > 0)
        S_r += -(t*H*MU/a0/NU)*(1 - 2*Y/(3*K*std::abs(S_r_c)));
    
    return S_r;
}

void next_layer_V(array<double, N_x> &V, array<double, N_x> &R, array<Point_x, N_x-1> &points) {
    
    for(int i = 1; i < N_x-1; i++) {
        V[i] = next_V(V[i], R[i-1], R[i], R[i+1], points[i-1], points[i]);
    }
    
    double P_l = 0;
    double P_r = P_right;
    
    V[0] = next_V_l(V[0], R[0], R[1], points[0], P_l);
    V[N_x-1] = next_V_r(V[N_x-1], R[N_x-2], R[N_x-1], points[N_x-2], P_r);
}

void next_Point_x(Point_x &point, double R_l, double R_r, double V_l, double V_r) {
    
    double e_r = get_e_R(V_l, V_r, R_l, R_r);
    double e_o = get_e_O(V_l, V_r, R_l, R_r);
    
    point.PO = next_PO(point, e_r, e_o);    
    point.PO = next_PO(point, e_r, e_o);
    point.SIGMA = -std::log(point.PO);
    point.SIGMA_R = point.SIGMA + point.S_R;
}

void next_layer(array<double, N_x> &V, array<double, N_x> &R, array<Point_x, N_x-1> &points) {
    array<double, N_x> R_new;
    
    next_layer_V(V, R, points);
    
    for(int i = 0; i < N_x; i++) {
        R_new[i] = next_R(R[i], V[i]);
    }
    
    for(int i = 0; i < N_x-1; i++) {
        double Rc_l = get_Rc(R_new[i], R[i]);
        double Rc_r = get_Rc(R_new[i+1], R[i+1]);
        
        next_Point_x(points[i], Rc_l, Rc_r, V[i], V[i+1]);
    }
    
    for(int i = 0; i < N_x; i++) R[i] = R_new[i];
}

void print_vals(array<double, N_x> &V, array<double, N_x> &R, array<Point_x, N_x-1> &points) {
    for(int i = 0; i < N_x-1; i++) {
        cout << R[i] << "\t";
        cout << V[i] << "\t";
        cout << points[i].PO << "\t";
        cout << points[i].SIGMA_R << "\t";
        cout << endl;
    }
    cout << R[N_x-1] << "\t";
    cout << V[N_x-1] << "\t";
    cout << points[N_x-2].PO << "\t";
    cout << P_right << "\t";
    cout << endl;
}


int main(int argc,char *argv[]) {
    
    double End = 0.001;
    array<double, N_x> R, V;
    array<Point_x, N_x-1> Points;
    
	if(argc > 1)
		End = atof(argv[1]);

    for(int i = 0; i < N_x; i++) {
        R[i] = (R1 + h*i)/H;
        V[i] = 0;
    }
    
    for(double Time = 0; Time < End; Time += t)
        next_layer(V, R, Points);
    
    print_vals(V, R, Points);
    
    return 0;
}


