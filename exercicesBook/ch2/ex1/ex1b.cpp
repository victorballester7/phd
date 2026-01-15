#include "../../libs/polylib.h"
#include <cmath>
#include <iostream>

using namespace std;

// map from [-1, 1] to [-1, 2]
double localToGlobal(double xi) { return 0.5 * (1.0 + 3 * xi); }



int main() {
  int *Q = new int[3];
  Q[0] = 4;
  Q[1] = 5;
  Q[2] = 6;

  double dxi_dx = 3./2.;


  double *w, *z;

  double sum = 0.0;
  for (int i = 0; i < 3; i++) {
    sum = 0.0;
    w = new double[Q[i]];
    z = new double[Q[i]];
    polylib::zwgll(z, w, Q[i]);

    for (int j = 0; j < Q[i]; j++)
      sum += w[j] * pow(localToGlobal(z[j]), 6);
    sum *= dxi_dx;

    cout << "Sum (" << Q[i] << ") = " << sum << endl;
    delete[] w;
    delete[] z;
  }
}
