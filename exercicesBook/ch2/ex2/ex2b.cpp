#include "../../libs/polylib.h"
#include <cmath>
#include <iostream>

using namespace std;

#define f(x) pow(x, 7)
#define df(x) 7 * pow(x, 6)

// map from [-1, 1] to [a, b]
double localToGlobal(double xi, double a, double b);

void dij(int Q, double *D, double *Dt);

double diff(int Q, double *D, int i, double a, double b);

int main() {
  double d = 0.0;

  double *D, *Dt, *w, *z;
  double a = 2.0;
  double b = 10.0;
  double dxi_dx = 2.0 / (b - a);

  for (int q = 7; q <= 9; q++) {
    D = new double[q * q];
    Dt = new double[q * q];
    dij(q, D, Dt);
    w = new double[q];
    z = new double[q];
    polylib::zwgll(z, w, q);

    for (int i = 0; i < q; i++) {
      d = diff(q, D, i, a, b);
      // plot with 16 digits of precision
      cout.precision(16);
      cout << "Sum (" << q << ", " << i << ") = " << d * dxi_dx  << ", " << df(localToGlobal(z[i], a, b)) << ", " << abs(d * dxi_dx - df(localToGlobal(z[i], a, b))) << endl;
    }
    delete[] D;
    delete[] Dt;
  }
}

double localToGlobal(double xi, double a, double b) {
  return 0.5 * (b + a + (b - a) * xi);
}

double globalToLocal(double x, double a, double b) {
  return 2 * (x - a) / (b - a) - 1;
}

void dij(int Q, double *D, double *Dt) {
  double *w, *z;
  w = new double[Q];
  z = new double[Q];
  polylib::zwgll(z, w, Q);
  polylib::Dgll(&D, &Dt, z, Q);

  delete[] w;
  delete[] z;
}

double diff(int Q, double *D, int i, double a, double b) {
  double sum = 0.0;
  double *w, *z;
  w = new double[Q];
  z = new double[Q];
  polylib::zwgll(z, w, Q);

  for (int j = 0; j < Q; j++) {
    sum += D[i * Q + j] * f(localToGlobal(z[j], a, b));
  }

  delete[] w;
  delete[] z;

  return sum;
}
