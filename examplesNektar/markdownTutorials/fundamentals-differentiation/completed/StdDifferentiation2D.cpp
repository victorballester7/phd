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
    cout << "========================================================" << endl;
    cout << "|   DIFFERENTIATION IN 2D ELEMENT in Standard Region   |" << endl;
    cout << "========================================================" << endl;
    cout << endl;
    cout << "Differentiate the function f(x1,x2) = (x1)^7*(x2)^9" << endl;
    cout << "in the standard quadrilateral element:" << endl;

    // Specify the number of quadrature points in both directions
    int nQuadPointsDir1 = 7;
    int nQuadPointsDir2 = 9;

    // Specify the type of quadrature points in both directions
    LibUtilities::PointsType quadPointsTypeDir1 =
        LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsType quadPointsTypeDir2 =
        LibUtilities::eGaussLobattoLegendre;

    // Declare variables (of type Array) to hold the quadrature zeros
    // and the values of the derivative at the quadrature points in both
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

    // Now you have the quadrature zeros and the differentiation matrix,
	// apply the Gaussian quadrature technique to differentiate the function
	// f(x_1,i,x_2,j) = x_1,i^7 * x2,j^9 on the standard
	// quadrilateral.  To do so, write a (double) loop which performs
// the summation.
//
// Store the solution in the matrices 'quadDerivsDir1' and 'quadDerivsDir2'

#if WITHSOLUTION
    for (size_t i = 0; i < nQuadPointsDir1; ++i)
    {
        for (size_t j = 0; j < nQuadPointsDir2; ++j)
        {
            for (size_t k = 0; k < nQuadPointsDir1; ++k)
            {
                quadDerivsDir1[i][j] += (*derivMatrixDir1)(i, k) *
                                        pow(quadZerosDir1[k], 7) *
                                        pow(quadZerosDir2[j], 9);
            }
            for (size_t k = 0; k < nQuadPointsDir2; ++k)
            {
                quadDerivsDir2[i][j] += (*derivMatrixDir2)(j, k) *
                                        pow(quadZerosDir1[i], 7) *
                                        pow(quadZerosDir2[k], 9);
            }
        }
    }
#endif

    // Compute the total error
    NekDouble error = 0.0;
    for (size_t i = 0; i < nQuadPointsDir1; ++i)
    {
        for (size_t j = 0; j < nQuadPointsDir2; ++j)
        {
            error +=
                fabs(quadDerivsDir1[i][j] -
                     7 * pow(quadZerosDir1[i], 6) * pow(quadZerosDir2[j], 9));
            error +=
                fabs(quadDerivsDir2[i][j] -
                     9 * pow(quadZerosDir1[i], 7) * pow(quadZerosDir2[j], 8));
        }
    }

    // Display the output
    cout << "\t q1 = " << nQuadPointsDir1 << ", q2 = " << nQuadPointsDir2;
    cout << ": Error = " << error << endl;
    cout << endl;
}
