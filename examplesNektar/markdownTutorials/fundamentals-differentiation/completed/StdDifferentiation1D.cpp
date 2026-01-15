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
    cout << "======================================================" << endl;
    cout << "|      DIFFERENTIATION IN A 1D STANDARD REGION       |" << endl;
    cout << "======================================================" << endl;

    cout << "Differentiate the function f(xi) = xi^7 in the " << endl;
    cout << "standard segment xi=[-1,1] using quadrature points" << endl;

    // Specify the number of quadrature points
    int nQuadPoints = 7;

    // Specify the type of quadrature points. This is done using the proper
    // Nektar++ syntax.
    LibUtilities::PointsType quadPointsType = LibUtilities::eGaussGaussLegendre;

    // Declare variables (of type Array) to hold the quadrature zeros
    // and the values of the derivative at the quadrature points
    Array<OneD, NekDouble> quadZeros(nQuadPoints);
    Array<OneD, NekDouble> quadDerivs(nQuadPoints);

    // Declare a pointer (to type NekMatrix<NekDouble>) to hold the
    // differentiation matrix
    DNekMatSharedPtr derivMatrix;

    // Calculate the quadrature zeros and the differentiation matrix.
    // This is done in 2 steps.

    // Step 1: Declare a PointsKey which uniquely defines the
    // quadrature points
    const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

    // Step 2: Using this key, the quadrature zeros and the differentiation
    // matrix can now be retrieved through the PointsManager in namespace
    // LibUtilities
    quadZeros   = (LibUtilities::PointsManager()[quadPointsKey])->GetZ();
    derivMatrix = (LibUtilities::PointsManager()[quadPointsKey])->GetD();

    // Now you have the quadrature zeros and the differentiation matrix,
	// apply the Gaussian quadrature technique to differentiate the function
	// f(xi) = xi^7 in the standard segment xi=[-1,1].  To do so, write a
	// loop which performs the summation for each quadrature point.
	//
	// In C++, a loop is implemented as:
	//     for(i = min; i < max; i++)
	//     {
	//         // do something
	//     }
	//
	// The function f(xi) can be evaluated using the command
	// 'pow(xi,7)' of the cmath standard library
	//
	// Store the solution in the array 'quadDerivs'

#if WITHSOLUTION
    for (size_t i = 0; i < nQuadPoints; ++i)
    {
        for (size_t j = 0; j < nQuadPoints; ++j)
        {
            quadDerivs[i] += (*derivMatrix)(i, j) * pow(quadZeros[j], 7);
        }
    }
#endif

    // Compute the total error
    NekDouble error = 0.0;
    for (size_t i = 0; i < nQuadPoints; ++i)
    {
        error += fabs(quadDerivs[i] - 7 * pow(quadZeros[i], 6));
    }

    // Display the output
    cout << "\t Q = " << nQuadPoints << ": Error = " << error << endl;

    // Now evaluate the derivatives for a quadrature order of Q = Q_max
    // where Q_max is the number of quadrature points required for an
    // exact evaluation of the derivative (calculate this value
    // analytically).  Check that the error should then be zero (up to
    // numerical precision).
}
