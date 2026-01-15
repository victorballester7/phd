#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>

using namespace Nektar;
using namespace std;

#define WITHSOLUTION 1

int main(int, char **)
{
    cout << "=========================================================" << endl;
    cout << "|     DIFFERENTIATION IN 2D ELEMENT in Local Region     |" << endl;
    cout << "=========================================================" << endl;
    cout << endl;
    cout << "Differentiate the function f(x1,x2) = x1^7 * x2^9 " << endl;
    cout << "in a local quadrilateral element:" << endl;

    // Specify the number of quadrature points in both directions
    int nQuadPointsDir1 = 8;
    int nQuadPointsDir2 = 10;

    // Specify the type of quadrature points in both directions
    LibUtilities::PointsType quadPointsTypeDir1 =
        LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsType quadPointsTypeDir2 =
        LibUtilities::eGaussLobattoLegendre;

    // Declare variables (of type Array) to hold the quadrature zeros
    // and the values of the derivatives at the quadrature points in both
    // directions
    Array<OneD, NekDouble> quadZerosDir1(nQuadPointsDir1);
    Array<OneD, NekDouble> quadZerosDir2(nQuadPointsDir2);
    Array<TwoD, NekDouble> quadDerivsDir1(nQuadPointsDir1, nQuadPointsDir2);
    Array<TwoD, NekDouble> quadDerivsDir2(nQuadPointsDir1, nQuadPointsDir2);

    // Declare pointers (to type NekMatrix<NekDouble>) to hold the
    // differentiation matrices
    DNekMatSharedPtr derivMatrixDir1;
    DNekMatSharedPtr derivMatrixDir2;

    // Calculate the GLL-quadrature zeros and the differentiation
    // matrices in both directions. This is done in 2 steps.

    // Step 1: Declare the PointsKeys which uniquely defines the
    // quadrature points
    const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1,
                                                    quadPointsTypeDir1);
    const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2,
                                                    quadPointsTypeDir2);

    // Step 2: Using this key, the quadrature zeros and differentiation
    // matrices can now be retrieved through the PointsManager
    quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
    derivMatrixDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetD();

    quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
    derivMatrixDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetD();

    // The local (straight-sided) quadrilateral element has the
    // following vertices:
    //
    // - Vertex A: (x1_A,x2_A) = (0,-1)
    // - Vertex B: (x1_A,x2_A) = (1,-1)
    // - Vertex C: (x1_A,x2_A) = (1,1)
    // - Vertex D: (x1_A,x2_A) = (0,0)
    //
    NekDouble x1_A = 0.0;
    NekDouble x2_A = -1.0;

    NekDouble x1_B = 1.0;
    NekDouble x2_B = -1.0;

    NekDouble x1_C = 1.0;
    NekDouble x2_C = 1.0;

    NekDouble x1_D = 0.0;
    NekDouble x2_D = 0.0;

    // Differentiate the function f(x1,x2) = x1^7 * x2^9 in a local
    // quadrilateral element (defined above). Use Gauss-Lobatto-Legendre
    // quadrature in both directions.
    //
    // Your code can be based on the previous exercise. However,
    // as we are calculating the derivatives of a function defined in
    // a local element rather than in a reference element, we have
    // to take into account the geometry of the element.
    //
    // Therefore, the implementation should be altered in two ways:
    //
    // (1) The quadrature zeros should be transformed to local
    //     coordinates to evaluate the function f(x1,x2)
    //
    // (2) Take into account the Jacobian matrix of the transformation
    //     between local and reference coordinates when evaluating
    //     the derivative. (Evaluate the expression for the Jacobian
    //     matrix analytically rather than using numerical
    //     differentiation and invert it)
    //
    // Store the solution in the matrices 'quadDerivsDir1' and 'quadDerivsDir2'

    // Apply the Gaussian quadrature technique to differentiate the
    // function f(x1,x2) = x1^7 * x2^9 in the standard
    // quadrilateral.  To do so, edit the (double) loop which performs
    // the summation.
    NekDouble x1_master, x2_master, x1_slave, x2_slave;
    NekDouble dx1dxi1, dx1dxi2, dx2dxi1, dx2dxi2;
    NekDouble dxi1dx1, dxi1dx2, dxi2dx1, dxi2dx2;
    NekDouble jacobian;
    NekDouble error = 0;

    for (size_t i = 0; i < nQuadPointsDir1; ++i)
    {
        for (size_t j = 0; j < nQuadPointsDir2; ++j)
        {
            // Compute the local coordinates of the quadrature point
            x1_master =
                x1_A * 0.25 * (1 - quadZerosDir1[i]) * (1 - quadZerosDir2[j]) +
                x1_B * 0.25 * (1 + quadZerosDir1[i]) * (1 - quadZerosDir2[j]) +
                x1_D * 0.25 * (1 - quadZerosDir1[i]) * (1 + quadZerosDir2[j]) +
                x1_C * 0.25 * (1 + quadZerosDir1[i]) * (1 + quadZerosDir2[j]);
            x2_master =
                x2_A * 0.25 * (1 - quadZerosDir1[i]) * (1 - quadZerosDir2[j]) +
                x2_B * 0.25 * (1 + quadZerosDir1[i]) * (1 - quadZerosDir2[j]) +
                x2_D * 0.25 * (1 - quadZerosDir1[i]) * (1 + quadZerosDir2[j]) +
                x2_C * 0.25 * (1 + quadZerosDir1[i]) * (1 + quadZerosDir2[j]);

            // Fill the Jacobian matrix analytically with the dx?dxi? terms
#if WITHSOLUTION
            dx1dxi1 = 0.25 * (1 - quadZerosDir2[j]) * (x1_B - x1_A) +
                      0.25 * (1 + quadZerosDir2[j]) * (x1_C - x1_D);
            dx2dxi2 = 0.25 * (1 - quadZerosDir1[i]) * (x2_D - x2_A) +
                      0.25 * (1 + quadZerosDir1[i]) * (x2_C - x2_B);
            dx1dxi2 = 0.25 * (1 - quadZerosDir1[i]) * (x1_D - x1_A) +
                      0.25 * (1 + quadZerosDir1[i]) * (x1_C - x1_B);
            dx2dxi1 = 0.25 * (1 - quadZerosDir2[j]) * (x2_B - x2_A) +
                      0.25 * (1 + quadZerosDir2[j]) * (x2_C - x2_D);
#endif

            // Compute the Jacobian determinant
            jacobian = dx1dxi1 * dx2dxi2 - dx1dxi2 * dx2dxi1;

            // Invert the Jacobian matrix to obtain the dxi?dx? terms
#if WITHSOLUTION
            dxi1dx1 = dx2dxi2 / jacobian;
            dxi2dx2 = dx1dxi1 / jacobian;
            dxi1dx2 = -dx1dxi2 / jacobian;
            dxi2dx1 = -dx2dxi1 / jacobian;
#endif

            for (size_t k = 0; k < nQuadPointsDir1; ++k)
            {
                // Compute the local coordinates of the quadrature point
                x1_slave = x1_A * 0.25 * (1 - quadZerosDir1[k]) *
                               (1 - quadZerosDir2[j]) +
                           x1_B * 0.25 * (1 + quadZerosDir1[k]) *
                               (1 - quadZerosDir2[j]) +
                           x1_D * 0.25 * (1 - quadZerosDir1[k]) *
                               (1 + quadZerosDir2[j]) +
                           x1_C * 0.25 * (1 + quadZerosDir1[k]) *
                               (1 + quadZerosDir2[j]);

                x2_slave = x2_A * 0.25 * (1 - quadZerosDir1[k]) *
                               (1 - quadZerosDir2[j]) +
                           x2_B * 0.25 * (1 + quadZerosDir1[k]) *
                               (1 - quadZerosDir2[j]) +
                           x2_D * 0.25 * (1 - quadZerosDir1[k]) *
                               (1 + quadZerosDir2[j]) +
                           x2_C * 0.25 * (1 + quadZerosDir1[k]) *
                               (1 + quadZerosDir2[j]);

                // Add the contribution of this quadrature point to
                // quadDerivsDir1[i][j] and quadDerivsDir2[i][j]
#if WITHSOLUTION
                quadDerivsDir1[i][j] += (*derivMatrixDir1)(i, k) *
                                        pow(x1_slave, 7) * pow(x2_slave, 9) *
                                        dxi1dx1;
                quadDerivsDir2[i][j] += (*derivMatrixDir1)(i, k) *
                                        pow(x1_slave, 7) * pow(x2_slave, 9) *
                                        dxi1dx2;
#endif
            }
            for (size_t k = 0; k < nQuadPointsDir2; ++k)
            {
                // Compute the local coordinates of the quadrature point
                x1_slave = x1_A * 0.25 * (1 - quadZerosDir1[i]) *
                               (1 - quadZerosDir2[k]) +
                           x1_B * 0.25 * (1 + quadZerosDir1[i]) *
                               (1 - quadZerosDir2[k]) +
                           x1_D * 0.25 * (1 - quadZerosDir1[i]) *
                               (1 + quadZerosDir2[k]) +
                           x1_C * 0.25 * (1 + quadZerosDir1[i]) *
                               (1 + quadZerosDir2[k]);

                x2_slave = x2_A * 0.25 * (1 - quadZerosDir1[i]) *
                               (1 - quadZerosDir2[k]) +
                           x2_B * 0.25 * (1 + quadZerosDir1[i]) *
                               (1 - quadZerosDir2[k]) +
                           x2_D * 0.25 * (1 - quadZerosDir1[i]) *
                               (1 + quadZerosDir2[k]) +
                           x2_C * 0.25 * (1 + quadZerosDir1[i]) *
                               (1 + quadZerosDir2[k]);

                // Add the contribution of this quadrature point to
                // quadDerivsDir1[i][j] and quadDerivsDir2[i][j]
#if WITHSOLUTION
                quadDerivsDir2[i][j] += (*derivMatrixDir2)(j, k) *
                                        pow(x1_slave, 7) * pow(x2_slave, 9) *
                                        dxi2dx2;
                quadDerivsDir1[i][j] += (*derivMatrixDir2)(j, k) *
                                        pow(x1_slave, 7) * pow(x2_slave, 9) *
                                        dxi2dx1;
#endif
            }

            error += fabs(quadDerivsDir1[i][j] -
                          7 * pow(x1_master, 6) * pow(x2_master, 9));
            error += fabs(quadDerivsDir2[i][j] -
                          9 * pow(x1_master, 7) * pow(x2_master, 8));
        }
    }

    // Display the averag error
    cout << "\t q1 = " << nQuadPointsDir1 << ", q2 = " << nQuadPointsDir2;
    cout << ": Average Error = ";
    cout << error / (nQuadPointsDir1 * nQuadPointsDir2) << endl;
    cout << endl;
}
