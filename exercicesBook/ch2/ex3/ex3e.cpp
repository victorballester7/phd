#include "../../libs/polylib.h"
#include <cmath>
#include <cstring>
#include <iostream>

#define f(x) pow(x, 7)

using namespace std;

extern "C" {
extern void dgetrf_(int *, int *, double(*), int *, int[], int *);
}
extern "C" {
extern void dgetrs_(char *, int *, int *, double(*), int *, int[], double(*),
                    int *, int *);
}

double localToGlobal(double xi, double a, double b);
void func(int Q, double *z, double *f, double a, double b);
double gaussIntegration(int Q, double *z, double *w, double *f1, double *f2,
                        double a, double b);
void buildPsiMatrix(int Q, int P, double *psiMatrix, double *z, double *w);
void buildPhiMatrix(int Q, int P, double *phiMatrix, double *z, double *w);
void buildM(int Q, int P, double *M, double *polyMatrix, double *z, double *w,
            double a, double b);
void buildf(int Q, int P, double *f, double *fquadPoints, double *bc,
            double *polyMatrix, double *z, double *w, double a, double b);
void solveSystem(int Q, int P, double *M, double *f);
void getUapprox(int Q, int P, double *u_approx, double *f, double *bc,
                double *polyMatrix);

int main() {
  int P = 8;
  int Q = 9;
  double a = 2.0;
  double b = 5.0;
  double *polyMatrix; // matrix of polynomials psi_p(z_i), for p = 0, 1, ..., P
                      // and i = 0, 1, ..., Q. psiMatrix[p*Q + i] = psi_p(z_i)
  double *u_approx;   // vector containung sum_{p=0}^{P} u_p psi_p(z_i) whose
                      // coordinates are the values of the function at the
                      // quadrature points
  double *M, *f, *fquadPoints;
  double *bc = new double[2];
  bc[0] = f(a);
  bc[1] = f(b);

  polyMatrix = new double[(P + 1) * Q];
  u_approx = new double[Q];
  M = new double[(P - 1) * (P - 1)];
  f = new double[P - 1];
  fquadPoints = new double[Q];

  // creation of quadrature points and weights
  double *z, *w;
  z = new double[Q];
  w = new double[Q];
  polylib::zwgll(z, w, Q);

  buildPhiMatrix(Q, P, polyMatrix, z, w);
  buildM(Q, P, M, polyMatrix, z, w, a, b);
  buildf(Q, P, f, fquadPoints, bc, polyMatrix, z, w, a, b);
  solveSystem(Q, P, M, f);
  getUapprox(Q, P, u_approx, f, bc, polyMatrix);

  for (int i = 0; i < Q; i++) {
    cout << "f(xi_" << i << ") = " << fquadPoints[i] << "\t" << u_approx[i]
         << "\t" << abs(fquadPoints[i] - u_approx[i]) << endl;
  }

  delete[] polyMatrix;
  delete[] M;
  delete[] f;
  delete[] fquadPoints;
  delete[] u_approx;
  delete[] z;
  delete[] w;
  return 0;
}

double localToGlobal(double xi, double a, double b) {
  return 0.5 * (b + a + (b - a) * xi);
}

void func(int Q, double *z, double *f, double a, double b) {
  for (int i = 0; i < Q; i++)
    f[i] = f(localToGlobal(z[i], a, b));
}

double gaussIntegration(int Q, double *z, double *w, double *f1, double *f2,
                        double a, double b) {
  double dx_dxi = 2. / (b - a);
  double integral = 0;
  for (int i = 0; i < Q; i++) {
    integral += w[i] * f1[i] * f2[i];
  }
  return integral * dx_dxi;
}

void buildPsiMatrix(int Q, int P, double *psiMatrix, double *z, double *w) {
  double *pjacob;
  pjacob = new double[Q];

  for (int p = 0; p < P + 1; p++) {
    if (p == 0) {
      for (int i = 0; i < Q; i++)
        psiMatrix[p * Q + i] = 0.5 * (1 - z[i]);
    } else if (p == P) {
      for (int i = 0; i < Q; i++)
        psiMatrix[p * Q + i] = 0.5 * (1 + z[i]);
    } else {
      polylib::jacobf(Q, z, pjacob, p - 1, 1, 1);
      for (int i = 0; i < Q; i++)
        psiMatrix[p * Q + i] = 0.25 * (1 - z[i]) * (1 + z[i]) * pjacob[i];
    }
  }
  delete[] pjacob;
}

void buildPhiMatrix(int Q, int P, double *phiMatrix, double *z, double *w) {
  double *pLegendre, *pLegendreDeriv;
  pLegendre = new double[Q];
  pLegendreDeriv = new double[Q];
  polylib::jacobf(Q, z, pLegendre, P, 0, 0);

  // compute derivative of Legendre polynomial
  polylib::jacobf(Q, z, pLegendreDeriv, P - 1, 1, 1);
  for (int i = 0; i < Q; i++) {
    pLegendreDeriv[i] *= 0.5 * (P + 1);
  }

  for (int p = 0; p < P + 1; p++) {
    for (int i = 0; i < Q; i++) {
      if (z[p] == z[i]) {
        phiMatrix[p * Q + i] = 1;
      } else {
        phiMatrix[p * Q + i] = (z[i] - 1) * (z[i] + 1) * pLegendreDeriv[i] /
                               (P * (P + 1) * pLegendre[p] * (z[p] - z[i]));
      }
    }
  }
  delete[] pLegendre;
  delete[] pLegendreDeriv;
}

void buildM(int Q, int P, double *M, double *polyMatrix, double *z, double *w,
            double a, double b) {
  for (int p = 0; p < P - 1; p++) {
    for (int q = 0; q < P - 1; q++)
      M[p * (P - 1) + q] = gaussIntegration(Q, z, w, polyMatrix + (p + 1) * Q,
                                            polyMatrix + (q + 1) * Q, a, b);
  }
}

void buildf(int Q, int P, double *f, double *fquadPoints, double *bc,
            double *polyMatrix, double *z, double *w, double a, double b) {
  // compute the values of the function f at the quadrature points
  func(Q, z, fquadPoints, a, b);

  double *f_aux = new double[Q];
  
  memcpy(f_aux, fquadPoints, Q * sizeof(double));

  for (int i = 0; i < Q; i++) {
    f_aux[i] -= bc[0] * polyMatrix[i];
    f_aux[i] -= bc[1] * polyMatrix[P * Q + i];
  }

  for (int p = 0; p < P - 1; p++)
    f[p] =
        gaussIntegration(Q, z, w, f_aux, polyMatrix + (p + 1) * Q, a, b);

  delete[] f_aux;
}

void solveSystem(int Q, int P, double *M, double *f) {
  int info;
  int *ipiv;
  int dim = P - 1;
  ipiv = new int[dim];
  char TRANS = 'T';
  int NRHS = 1;

  dgetrf_(&dim, &dim, M, &dim, ipiv, &info);
  dgetrs_(&TRANS, &dim, &NRHS, M, &dim, ipiv, f, &dim, &info);

  delete[] ipiv;
}

void getUapprox(int Q, int P, double *u_approx, double *f, double *bc,
                double *polyMatrix) {
  for (int i = 0; i < Q; i++) {
    u_approx[i] = 0;
    for (int p = 0; p < P - 1; p++) {
      u_approx[i] += f[p] * polyMatrix[(p + 1) * Q + i];
    }
    u_approx[i] += bc[0] * polyMatrix[i];
    u_approx[i] += bc[1] * polyMatrix[P * Q + i];
  }
}
