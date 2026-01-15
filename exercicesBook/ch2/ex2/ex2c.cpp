#include "../../libs/polylib.h"
#include <cmath>
#include <iostream>

using namespace std;

#define f(x) cos(x)

// map from [-1, 1] to [a, b]
double localToGlobal(double xi, double a, double b);

void dij(int Q, double *D, double *Dt);

int main() {
  double d = 0.0;

  double *D, *Dt, *w, *z;
  double a = 0.0;
  double b = M_PI_2;

  cout.precision(16);
  for (int q = 2; q <= 8; q++) {
    D = new double[q * q];
    Dt = new double[q * q];
    dij(q, D, Dt);
    w = new double[q];
    z = new double[q];
    polylib::zwgll(z, w, q);

    d = 0.0;
    for (int i = 0; i < q; i++) {
      for (int j = 0; j < q; j++) {
         d += w[i] * D[i * q + j] * f(localToGlobal(z[j], a, b));
      }
    }
    cout << "Sum (" << q << ") = " << -d << ", " << abs(-d - 1) << endl;

    delete[] D;
    delete[] Dt;
  }
}

double localToGlobal(double xi, double a, double b) {
  return 0.5 * (b + a + (b - a) * xi);
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

