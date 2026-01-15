#include <cmath>
#include "../../libs/polylib.h"
#include <iostream>

using namespace std;


int main(){
  int* Q = new int[3];
  Q[0] = 4;
  Q[1] = 5;
  Q[2] = 6;

  double* w, *z;
  double sum = 0.0;
  for(int i = 0; i < 3; i++){
    sum = 0.0;
    w = new double[Q[i]];
    z = new double[Q[i]];
    polylib::zwgll(z, w, Q[i]);

    for(int j = 0; j < Q[i]; j++)
      sum += w[j] * pow(z[j], 6);

    cout << "Sum (" << Q[i] << ") = " << sum << endl;
    delete[] w;
    delete[] z;
  }

}
