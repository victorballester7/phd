#include "../../libs/polylib.h"
#include <cmath>
#include <cstring>
#include <iostream>

#define f(x) sin(x)
#define map(e, p) mapArray[(e) * (P + 1) + (p)]

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
void buildMapArray(int Nel, int P, int *map);
void buildMe(int Q, int P, double *M, double *polyMatrix, double *z, double *w,
             double a, double b);
void buildMg(int Nel, int Q, int P, int dim, double *Mg, int *mapArray,
             double *polyMatrix, double *z, double *w, double A, double B);
void buildfe(int Nel, int e, int Q, int P, double *f, double *fquadPoints, double *bc,
             double *polyMatrix, double *z, double *w, double a, double b);
void buildfg(int Nel, int Q, int P, int dim, double *fg, int *mapArray,
             double *fquadPoints, double *bc, double *polyMatrix, double *z,
             double *w, double A, double B);
void solveSystem(int dim, double *Mg, double *fg);
void getUapprox(int Nel, int Q, int P, double *u_approx, double *fg, double *bc,
                int *mapArray, double *polyMatrix);

int main() {
  int Nel = 10;
  int P = 8;
  int Q = 10;
  int dim = Nel * P - 1;
  double A = 0.0;
  double B = 10.0;
  double *polyMatrix; // matrix of polynomials psi_p(z_i), for p = 0, 1, ..., P
                      // and i = 0, 1, ..., Q. psiMatrix[p*Q + i] = psi_p(z_i)
  double *u_approx;   // vector containung sum_{p=0}^{P} u_p psi_p(z_i) whose
                      // coordinates are the values of the function at the
                      // quadrature points
  double *Mg, *fg, *fquadPoints;
  int *mapArray;
  double *bc = new double[2];
  bc[0] = f(A);
  bc[1] = f(B);

  polyMatrix = new double[(P + 1) * Q];
  mapArray = new int[Nel * (P + 1)];
  u_approx = new double[Nel * (Q - 1) + 1];
  Mg = new double[dim * dim];
  fg = new double[dim];
  fquadPoints = new double[Q];

  // creation of quadrature points and weights
  double *z, *w;
  z = new double[Q];
  w = new double[Q];
  polylib::zwgll(z, w, Q);

  buildPsiMatrix(Q, P, polyMatrix, z, w);
  buildMapArray(Nel, P, mapArray);
  buildMg(Nel, Q, P, dim, Mg, mapArray, polyMatrix, z, w, A, B);
  // print non-zero elements of Mg
  // for (int i = 0; i < dim; i++) {
  //   for (int j = 0; j < dim; j++) {
  //     if (abs(Mg[i * dim + j]) > 1e-10)
  //       cout << Mg[i * dim + j] << endl;
  //   }
  // }
  buildfg(Nel, Q, P, dim, fg, mapArray, fquadPoints, bc, polyMatrix, z, w, A,
          B);
  // print vector
  // for (int i = 0; i < dim; i++) {
  //   cout << fg[i] << endl;
  // }
  solveSystem(dim, Mg, fg);
  getUapprox(Nel, Q, P, u_approx, fg, bc, mapArray, polyMatrix);

  for (int e = 0; e < Nel; e++) {
    for (int i = 0; i < Q; i++) {
      if (i == Q - 1 && e != Nel - 1)
        continue;
      func(Q, z, fquadPoints, A + e * (B - A) / Nel,
           A + (e + 1) * (B - A) / Nel);
      cout << "f(xi_" << e << "," << i << ") = " << fquadPoints[i] << "\t"
           << u_approx[e * (Q - 1) + i] << "\t"
           << abs(fquadPoints[i] - u_approx[e * (Q - 1) + i]) << endl;
    }
  }

  delete[] polyMatrix;
  delete[] Mg;
  delete[] fg;
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
  double dx_dxi = (b - a) / 2.0;
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

void buildMapArray(int Nel, int P, int *map) {
  for (int e = 0; e < Nel; e++) {
    for (int p = 0; p < P + 1; p++) {
      if (e == 0 && p == 0)
        map[e * (P + 1) + p] = Nel * P - 1;
      else if (e == Nel - 1 && p == P)
        map[e * (P + 1) + p] = Nel * P;
      else
        map[e * (P + 1) + p] = e * P + p - 1;
    }
  }
}

void buildMe(int Q, int P, double *M, double *polyMatrix, double *z, double *w,
             double a, double b) {
  for (int p = 0; p < P + 1; p++) {
    for (int q = 0; q < P + 1; q++)
      M[p * (P + 1) + q] = gaussIntegration(Q, z, w, polyMatrix + p * Q,
                                            polyMatrix + q * Q, a, b);
  }
}

void buildMg(int Nel, int Q, int P, int dim, double *Mg, int *mapArray,
             double *polyMatrix, double *z, double *w, double A, double B) {
  double *Me;
  Me = new double[(P + 1) * (P + 1)];
  double el_length = (B - A) / Nel;

  for (int i = 0; i < dim * dim; i++) {
    Mg[i] = 0;
  }

  for (int e = 0; e < Nel; e++) {
    buildMe(Q, P, Me, polyMatrix, z, w, A + e * el_length,
            A + (e + 1) * el_length);
    for (int p = 0; p < P + 1; p++) {
      for (int q = 0; q < P + 1; q++) {
        if (e == 0 && (p == 0 || q == 0))
          continue;
        if (e == Nel - 1 && (p == P || q == P))
          continue;
        Mg[map(e, p) * dim + map(e, q)] += Me[p * (P + 1) + q];
      }
    }
  }
  delete[] Me;
}

void buildfe(int Nel, int e, int Q, int P, double *f, double *fquadPoints, double *bc,
             double *polyMatrix, double *z, double *w, double a, double b) {
  // compute the values of the function f at the quadrature points
  func(Q, z, fquadPoints, a, b);

  double *f_aux = new double[Q];

  memcpy(f_aux, fquadPoints, Q * sizeof(double));

  if (e == 0) {
    for (int i = 0; i < Q; i++) {
      f_aux[i] -= bc[0] * polyMatrix[i];
    }
  } else if (e == Nel - 1) {
    for (int i = 0; i < Q; i++) {
      f_aux[i] -= bc[1] * polyMatrix[P * Q + i];
    }
  }

  for (int p = 0; p < P + 1; p++)
    f[p] = gaussIntegration(Q, z, w, f_aux, polyMatrix + p * Q, a, b);

  delete[] f_aux;
}

void buildfg(int Nel, int Q, int P, int dim, double *fg, int *mapArray,
             double *fquadPoints, double *bc, double *polyMatrix, double *z,
             double *w, double A, double B) {
  double *fe;
  fe = new double[P + 1];
  double el_length = (B - A) / Nel;

  for (int i = 0; i < dim; i++) {
    fg[i] = 0;
  }

  for (int e = 0; e < Nel; e++) {
    buildfe(Nel, e, Q, P, fe, fquadPoints, bc, polyMatrix, z, w, A + e * el_length,
            A + (e + 1) * el_length);
    if (e == 0) {
      for (int p = 1; p < P + 1; p++) {
        fg[map(e, p)] += fe[p];
      }
    } else if (e == Nel - 1) {
      for (int p = 0; p < P; p++) {
        fg[map(e, p)] += fe[p];
      }
    } else {
      for (int p = 0; p < P + 1; p++) {
        fg[map(e, p)] += fe[p];
      }
    }
  }
  delete[] fe;
}

void solveSystem(int dim, double *M, double *f) {
  int info;
  int *ipiv;
  ipiv = new int[dim];
  char TRANS = 'T';
  int NRHS = 1;

  dgetrf_(&dim, &dim, M, &dim, ipiv, &info);
  dgetrs_(&TRANS, &dim, &NRHS, M, &dim, ipiv, f, &dim, &info);

  delete[] ipiv;
}

void getUapprox(int Nel, int Q, int P, double *u_approx, double *fg, double *bc,
                int *mapArray, double *polyMatrix) {
  for (int e = 0; e < Nel; e++) {
    for (int i = 0; i < Q; i++) {
      if (i == Q - 1 && e != Nel - 1)
        continue;
      u_approx[e * (Q - 1) + i] = 0;
      if (e == 0) {
        for (int p = 1; p < P + 1; p++) {
          u_approx[e * (Q - 1) + i] += fg[map(e, p)] * polyMatrix[p * Q + i];
        }
        u_approx[e * (Q - 1) + i] += bc[0] * polyMatrix[i];
      } else if (e == Nel - 1) {
        for (int p = 0; p < P; p++) {
          u_approx[e * (Q - 1) + i] += fg[map(e, p)] * polyMatrix[p * Q + i];
        }
        u_approx[e * (Q - 1) + i] += bc[1] * polyMatrix[P * Q + i];
      } else {
        for (int p = 0; p < P + 1; p++) {
          u_approx[e * (Q - 1) + i] += fg[map(e, p)] * polyMatrix[p * Q + i];
        }
      }
    }
  }
}
