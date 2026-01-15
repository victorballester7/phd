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
    cout << "|      INTEGRATION ON 2D ELEMENT in Local Region          |"<< endl;
    cout << "==========================================================="<< endl;
    cout << endl;
    cout << "Integrate the function f(x1,x2) = x1^12 * x2^14 "   << endl; 
    cout << "on a local quadrilateral element:" << endl;

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
    
    
    // The local (straight-sided) quadrilateral element has the
    // following vertices:
    //
    // - Vertex A: (x1_A,x2_A) = (0,-1)
    // - Vertex B: (x1_A,x2_A) = (1,-1)
    // - Vertex C: (x1_A,x2_A) = (1,1)
    // - Vertex D: (x1_A,x2_A) = (0,0)
    //
    NekDouble x1_A =  0.0;
    NekDouble x2_A = -1.0;

    NekDouble x1_B =  1.0;
    NekDouble x2_B = -1.0;

    NekDouble x1_C =  1.0;
    NekDouble x2_C =  1.0;

    NekDouble x1_D =  0.0;
    NekDouble x2_D =  0.0;       
    
    // Integrate the function f(x1,x2) = x1^12 * x2^14 on a local
    // quadrilateral element (defined above). Use 7th order
    // Gauss-Lobatto-Legendre quadrature in both direction.
    //
    // Your code can be based on the previous exercise. However,
    // as we are calculating the integral of a function defined on
    // a local element rather than on a reference element, we have
    // to take into account the geometry of the element.
    //
    // Therefore, the implementation should be altered in two ways:
    //
    // (1) The quadrature zeros should be transformed to local
    //     coordinates to evaluate the integrand f(x1,x2)
    //
    // (2) Take into account the Jacobian of the transformation
    //     between local and reference coordinates when evaluating
    //     the integral. (Evaluate the expression for the Jacobian
    //     analytically rather than using numerical
    //     differentiation)
    //

    // Store the solution in the variable 'result'
    NekDouble result = 0.0;

    // Apply the Gaussian quadrature technique to integrate the
    // function f(x1,x2) = x1^12 * x2^14 on the standard
    // quadrilateral.  To do so, edit the (double) loop which performs
    // the summation.
    double x1;
    double x2;
    double jacobian;
    
    for(int i = 0; i < nQuadPointsDir1; i++)
    { 
        for(int j = 0; j < nQuadPointsDir2; j++)
        {
            // Calculate the local coordinates of the quadrature zeros
            // using the mapping from reference to local element
            x1 = x1_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                 x1_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                 x1_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                 x1_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 
            
            x2 = x2_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                 x2_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                 x2_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                 x2_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 
            
#if WITHSOLUTION
            // Analytically evaluate the Jacobian
            jacobian = (0.25 * (1-quadZerosDir2[j]) * (x1_B-x1_A)  + 
                        0.25 * (1+quadZerosDir2[j]) * (x1_C-x1_D)) *
                       (0.25 * (1-quadZerosDir1[i]) * (x2_D-x2_A)  +
                        0.25 * (1+quadZerosDir1[i]) * (x2_C-x2_B)) - 
                       (0.25 * (1-quadZerosDir1[i]) * (x1_D-x1_A)  +
                        0.25 * (1+quadZerosDir1[i]) * (x1_C-x1_B)) *
                       (0.25 * (1-quadZerosDir2[j]) * (x2_B-x2_A)  +
                        0.25 * (1+quadZerosDir2[j]) * (x2_C-x2_D));
#endif
            
            result += pow(x1,12) *  pow(x2,14) * 
                quadWeightsDir1[i] * quadWeightsDir2[j] * fabs(jacobian);
        }
    }
    
    // Display the output
    NekDouble exactResult = 1.0/195.0+1.0/420.0;
    cout << "\t Error = "<< fabs(result-exactResult) << endl;
}

