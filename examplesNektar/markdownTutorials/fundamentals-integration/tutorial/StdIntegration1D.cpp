#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Polylib/Polylib.h>

using namespace Nektar;
using namespace std;

int main(int, char **) {
  cout << "======================================================" << endl;
  cout << "|        INTEGRATION ON A 1D STANDARD REGION         |" << endl;
  cout << "======================================================" << endl;

  cout << "Integrate the function f(xi) = xi^12 on the standard " << endl;
  cout << "segment xi=[-1,1] with Gaussian quadrature" << endl;

  // Specify the number of quadrature points
  int nQuadPoints = 4;

  // Specify the type of quadrature points. This is done using the proper
  // Nektar++ syntax.
  LibUtilities::PointsType quadPointsType = LibUtilities::eGaussGaussLegendre;

  // Declare variables (of type Array) to hold the quadrature zeros
  // and weights
  Array<OneD, NekDouble> quadZeros(nQuadPoints);
  Array<OneD, NekDouble> quadWeights(nQuadPoints);

  // Calculate the quadrature zeros and weights. This is done in 2
  // steps.  Step 1: Declare a PointsKey which uniquely defines the
  // quadrature points
  const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

  // Step 2: Using this key, the quadrature zeros and weights can
  // now be retrieved through the PointsManager in namespace LibUtilities
  quadZeros = (LibUtilities::PointsManager()[quadPointsKey])->GetZ();
  quadWeights = (LibUtilities::PointsManager()[quadPointsKey])->GetW();

  // Now you have the quadrature zeros and weight, apply the
  // Gaussian quadrature technique to integrate the function f(xi) =
  // xi^12 on the standard segment xi=[-1,1].  To do so, write a
  // loop which performs the summation.
  //
  // In C++, a loop is implemented as:
  //     for(i = min; i < max; i++)
  //     {
  //         // do something
  //     }
  //
  // The function f(xi) can be evaluated using the command
  // 'pow(xi,12)' of the cmath standard library
  //
  // Store the solution in the variable 'result'
  NekDouble result = 0.0;
  for (int i = 0; i < nQuadPoints; i++) {
    result += pow(quadZeros[i], 12) * quadWeights[i];
  }
  //                 ==> Write your code here <==

  // Display the output
  NekDouble exactResult = 2.0 / 13.0;
  cout << "\t Q = " << nQuadPoints << ": Error = " << fabs(result - exactResult)
       << endl;

  // Now evaluate the integral for a quadrature order of Q = Q_max
  // where Q_max is the number of quadrature points required for an
  // exact evaluation of the integral (calculate this value
  // analytically).  Check that the error should then be zero (up to
  // numerical precision).
}
