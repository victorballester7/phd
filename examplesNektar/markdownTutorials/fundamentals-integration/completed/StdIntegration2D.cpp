#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
using namespace Nektar;
using namespace std;

#define WITHSOLUTION 1

int main(int, char **)
{
    cout << "==========================================================="<< endl;
    cout << "|    INTEGRATION ON 2D ELEMENT in Standard Region         |"<< endl;
    cout << "==========================================================="<< endl;
    cout << endl;
    cout << "Integrate the function f(x1,x2) = (x1)^12*(x2)^14" << endl;
    cout << "on the standard quadrilateral element:" << endl;
    
    // Specify the number of quadrature points in both directions
    int nQuadPointsDir1 = 6;
    int nQuadPointsDir2 = 7;
    
    // Specify the type of quadrature points in both directions 
    LibUtilities::PointsType quadPointsTypeDir1 =
        LibUtilities::eGaussLobattoLegendre;
    LibUtilities::PointsType quadPointsTypeDir2 =
        LibUtilities::eGaussLobattoLegendre;
    
    // Declare variables (of type Array) to hold the quadrature
    // zeros and weights
    Array<OneD, NekDouble> quadZerosDir1  (nQuadPointsDir1);
    Array<OneD, NekDouble> quadWeightsDir1(nQuadPointsDir1);
    Array<OneD, NekDouble> quadZerosDir2  (nQuadPointsDir2);
    Array<OneD, NekDouble> quadWeightsDir2(nQuadPointsDir2);
    
    // Calculate the GLL-quadrature zeros and weights in both
    // directions. This is done in 2 steps.
    //
    // Step 1: Declare the PointsKeys which uniquely defines the
    // quadrature points
    const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1,
                                                    quadPointsTypeDir1);
    const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2,
                                                    quadPointsTypeDir2);
    
    // Step 2: Using this key, the quadrature zeros and weights
    // can now be retrieved through the PointsManager
    quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
    quadWeightsDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetW();
    
    quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
    quadWeightsDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetW();
    
    // Now you have the quadrature zeros and weight, apply the
    // Gaussian quadrature technique to integrate the function
    // f(x_1,i,x_2,j) = x_1,i^12 * x2,j^14 on the standard
    // quadrilateral.  To do so, write a (double) loop which performs
    // the summation.
    //
    // Store the solution in the variable 'result'
    NekDouble result = 0.0;
    
#if WITHSOLUTION
    for(int i = 0; i < nQuadPointsDir1; i++)
    { 
        for(int j = 0; j < nQuadPointsDir2; j++)
        {
            result += pow(quadZerosDir1[i],12)*pow(quadZerosDir2[j],14)*
                quadWeightsDir1[i]*quadWeightsDir2[j];
        }
    }
#endif

    // Display the output
    NekDouble exactResult = 4.0/195.0;
    cout << "\t q1 = " << nQuadPointsDir1 << ", q2 = " << nQuadPointsDir2 ;
    cout << ": Error = "<< fabs(result-exactResult) << endl;
    cout << endl;
}

