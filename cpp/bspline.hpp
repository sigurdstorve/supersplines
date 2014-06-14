#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <stdexcept>

// Workaround for MinGW bug
#pragma GCC optimize ("-fno-ipa-cp-clone")

namespace bspline_storve {

// Simple implementation of functionality like NumPy's linspace()
void linspace(float left, float right, int n, std::vector<float>& res) {
    if (n < 0) {
        return;
    }
    res.resize(n);
    if (n == 0) {
        return;
    } else if (n == 1) {
        res[0] = left;
    } else {
        for (int i = 0; i < n; i++) {
            res[i] = left+i*(right-left)/(n-1);
        }
    }
}

// Test if floating point number is zero.
template <typename T>
bool _floatIsZero(T value, double eps=1e-6) {
    return std::abs(value) < eps;
}

//Return num/dev with the special rule that 0/0 is 0.
template <typename T>
T _specialDiv(T num, T den) {
    if (_floatIsZero(num) && _floatIsZero(den)) {
        return 0.0;
    } else {
        return num / den;
    }
}

// Compute B-splines directly using recursive definition.        
//  j:  B-spline basis no.
//  p:  Polynomial degree
//  knots: knot vector (must support indexing with [])
float bsplineBasis(int j, int p, float x, const std::vector<float>& knots) {
    if (p == 0 ) {
        if ( (x >= knots[j]) && (x < knots[j+1]) ) {
            return 1.0f;
        } else {
            return 0.0f;
        }
    } else {
        float left = _specialDiv((x-knots[j])*bsplineBasis(j,p-1,x,knots), knots[j+p]-knots[j]);
        float right = _specialDiv((knots[j+1+p]-x)*bsplineBasis(j+1,p-1,x,knots), knots[j+1+p]-knots[j+1]);
        return left + right;
    }
}

/*
    Find the knot interval index (mu) for a given value x
    This is a prerequisite for alg 2.20 and 2.21.
    Optional to supply a suggestion to try first.
    Throws on error.
*/
int getMu(float x, const std::vector<float>& knots, int mu=-1) {
    if (mu != -1 && (mu <= knots.size()-2) && (mu >= 0)) {
        if (x >= knots[mu] && x < knots[mu+1]) {
            return mu;
        }
    }
    for (int mu = 0; mu < (knots.size()-1); mu++) {
        if (x >= knots[mu] && x < knots[mu+1]) {
            return mu;
        }
    }
    throw std::runtime_error("Illegal knot vector or x value");
}

// Compute B-spline matrix from theorem 2.18
Eigen::MatrixXf bsplineMatrix(int k, int mu, float x, const std::vector<float>& knots) {
    if (k <= 0) throw std::runtime_error("k must be greater than zero");
    Eigen::MatrixXf res(k, k+1);
    for (int row = 0; row <k; row++) {
        for (int col = 0; col < k+1; col++) {
            res(row,col) = 0.0f;
        }
    }    
    for (int row = 0; row < k; row++) {
        float common_denom = knots[mu+(row+1)] - knots[mu+(row+1)-k];
        res(row,row) = _specialDiv(knots[mu+(row+1)] - x, common_denom);
        res(row,row+1) = _specialDiv(x - knots[mu + (row+1) - k], common_denom);
    }
    return res;
}

// Compute B-spline vector using algorithm 2.19/2.20
// Returns (B_{\mu-p,p}(x),...,B_{\mu,p}(x))
Eigen::VectorXf bsplineVector(int p, int mu, float x, const std::vector<float>& knots) {
    Eigen::MatrixXf res = bsplineMatrix(1, mu, x, knots);
    for (int k = 2; k <= p; k++) {
        res *= bsplineMatrix(k, mu, x, knots);
    }
    return res;
}

// Direkte
float B0(int j, float x, const std::vector<float>& knots) {
    if (knots[j] <= x && x < knots[j+1]) {
        return 1.0f;
    } else {
        return 0.0f;
    }
}

float B1(int j, float x, const std::vector<float>& knots) {
    if (knots[j] <= x && x < knots[j+1]) {
        // *B_{j,0}
        return (
            (x-knots[j])/(knots[j+1]-knots[j])
        );
    } else if (knots[j+1] <= x && x < knots[j+2]) {
        // *B_{j+1,0}
        return (
            (knots[j+2]-x) / (knots[j+2]-knots[j+1])
        );
    } else {
        return 0.0f;
    } 
}

float B2(int j, float x, const std::vector<float>& knots) {
    if (knots[j] <= x && x < knots[j+1]) {
        // *B_{j,0}
        return (
            ((x-knots[j])/(knots[j+2]-knots[j]))*((x-knots[j])/(knots[j+1]-knots[j]))
        );
    } else if (knots[j+1] <= x && x < knots[j+2]) {
        // *B_{j+1,0}
        return (
            ((x-knots[j])/(knots[j+2]-knots[j]))*((knots[j+2]-x)/(knots[j+2]-knots[j+1]))
            + ((knots[j+3]-x)/(knots[j+3]-knots[j+1]))*((x-knots[j+1])/(knots[j+2]-knots[j+1]))
        );
    } else if (knots[j+2] <= x && x < knots[j+3]) {
        // *B_{j+2,0}
        return (
            ((knots[j+3]-x)/(knots[j+3]-knots[j+1]))*((knots[j+3]-x)/(knots[j+3]-knots[j+2]))
        );
    } else {
        return 0.0f;
    }
}

float B3(int j, float x, const std::vector<float>& knots) {
    if (knots[j] <= x && x < knots[j+1]) {
        // *B_{j,0}
        return (
            ((x-knots[j])/(knots[j+3]-knots[j]))*((x-knots[j])/(knots[j+2]-knots[j]))*((x-knots[j])/(knots[j+1]-knots[j]))
        );
    } else if (knots[j+1] <= x && x < knots[j+2]) {
        // *B_{j+1,0}
        return (
            ((x-knots[j])/(knots[j+3]-knots[j]))*((x-knots[j])/(knots[j+2]-knots[j]))*((knots[j+2]-x)/(knots[j+2]-knots[j+1]))
            + ((x-knots[j])/(knots[j+3]-knots[j]))*((knots[j+3]-x)/(knots[j+3]-knots[j+1]))*((x-knots[j+1])/(knots[j+2]-knots[j+1]))
            + ((knots[j+4]-x)/(knots[j+4]-knots[j+1]))*((x-knots[j+1])/(knots[j+3]-knots[j+1]))*((x-knots[j+1])/(knots[j+2]-knots[j+1]))
        );
    } else if (knots[j+2] <= x && x < knots[j+3]) {
        // *B_{j+2,0}
        return (
            ((x-knots[j])/(knots[j+3]-knots[j]))*((knots[j+3]-x)/(knots[j+3]-knots[j+1]))*((knots[j+3]-x)/(knots[j+3]-knots[j+2]))
            + ((knots[j+4]-x)/(knots[j+4]-knots[j+1]))*((x-knots[j+1])/(knots[j+3]-knots[j+1]))*((knots[j+3]-x)/(knots[j+3]-knots[j+2]))
            + ((knots[j+4]-x)/(knots[j+4]-knots[j+1]))*((knots[j+4]-x)/(knots[j+4]-knots[j+2]))*((x-knots[j+2])/(knots[j+3]-knots[j+2]))
        );
    } else if (knots[j+3] <= x && x < knots[j+4]) {
        // *B_{j+3,0}
        return (
            ((knots[j+4]-x)/(knots[j+4]-knots[j+1]))*((knots[j+4]-x)/(knots[j+4]-knots[j+2]))*((knots[j+4]-x)/(knots[j+4]-knots[j+3]))
        );
    } else {
        return 0.0f;
    }
}

/*
Create a p+1-regular uniform knot vector for
a given number of control points
Throws if n is too small
endHack: Make the last knot bigger than the one before it, which has the
effect of making it non-zero in last knot value.
*/
std::vector<float> uniformRegularKnotVector(int numPoints,
                                            int degree,
                                            float tStart=0.0f,
                                            float tEnd=1.0f,
                                            bool endHack=false) {

    // The minimum length of a degree+1-regular knot vector is 2*(degree+1)
    if (numPoints < degree+1) {
        throw std::runtime_error("Too small n for a uniform regular knot vector");
    }
    
    std::vector<float> knots;
    // degree+1 copies of tStart left and degree+1 copies of tEnd right
    // but one of each in linspace
    for (int i = 0; i < degree; i++) {
        knots.push_back(tStart);
    }
    std::vector<float> temp;
    linspace(tStart, tEnd, numPoints+1-degree, temp);
    for (float t : temp) {
        knots.push_back(t);
    }
    for (int i = 0; i < degree; i++) {
        knots.push_back(tEnd);
    }
    
    if (endHack) {
        int numKnots = knots.size();
        knots[numKnots-1] += 1.0f;
    }
    return knots;
}

/*
Return the control point abscissa for the control polygon
of a one-dimensional spline.
*/
std::vector<float> controlPolygon(int p, const std::vector<float>& knots) {
    int n = knots.size() - p - 1;
    std::vector<float> abscissas;
    for (int i = 0; i < n; i++) { 
        float sum = 0.0f;
        for (int j = 1; j <= p; j++) {
            sum += knots[i+j];
        }
        abscissas.push_back(sum / p);
    }
    return abscissas;
}

/*
Compute points on a spline function using the straightforward
implementation of the recurrence relation for the B-splines.
*/
std::vector<float> renderSpline(int p,
                                const std::vector<float>& knots,
                                const std::vector<float>& controlPoints,
                                const std::vector<float>& ts) {
    std::vector<float> ys;
    for (float t : ts) {
        float y = 0.0f; 
        for (int j = 0; j < controlPoints.size(); j++) {
            y += bsplineBasis(j, p, t, knots)*controlPoints[j];
        }
        ys.push_back(y);
    }
    return ys;
}

/*
Approximate given data points by a spline in a 
given spline space through least-squares.
*/
std::vector<float> leastSquaresFit(const std::vector<float>& xs,
                                   const std::vector<float>& ys,
                                   const std::vector<float>& knots,
                                   int p,
                                   const std::vector<float>& _weights = {}) {
    int numKnots = knots.size();
    int numDataPoints = xs.size();
    //assert(ys.size() == numDataPoints);
   
    // The number of knots determines the number of control points
    // to use in the approximation.
    int numApproxPoints = numKnots - p - 1;
   
    std::vector<float> weights(numDataPoints);
    if (_weights.size() == 0) {
        for (int i = 0; i < numDataPoints; i++) weights[i] = 1.0f;
    } else {
       weights = _weights;
    }
    
    // Create matrix A
    Eigen::MatrixXf A(numDataPoints, numApproxPoints);
    for (int row = 0; row < numDataPoints; row++) {
        for (int col = 0; col < numApproxPoints; col++) {
            A(row, col) = std::sqrt(weights[row])*bsplineBasis(col, p, xs[row], knots);
        }
    }
        
    // Create vector b
    Eigen::VectorXf b(numDataPoints);
    for (int i = 0; i < numDataPoints; i++) {
        b(i) = std::sqrt(weights[i])*ys[i];
    }
    
    // Compute least-squares coefficients
    Eigen::VectorXf temp = (A.transpose()*A).ldlt().solve(A.transpose()*b);
    assert(temp.size() == numApproxPoints);
    std::vector<float> res(numApproxPoints);
    for (int i = 0; i < numApproxPoints; i++) {
        res[i] = temp(i);
    }
    return res;
}

// Compute the control points of the derivative of a given spline.
// The derivatoive of a degree-p spline with n B-spline coefficients
// is a degree-(p-1) spline with n+1 control points and can be
// expressed on the same knot vector.
std::vector<float> computeDerivativeCoeffs(int degree,
                                           const std::vector<float>& coeffs,
                                           const std::vector<float>& knots) {
    if (degree <= 0) throw std::runtime_error("Cannot differentiate a spline of degree lower than 1");
    int n = coeffs.size();
    std::vector<float> derCoeffs(n+1);
    // Handle boundary cases
    for (int i = 0; i < n+1; i++) {
        float value;
        if (i == 0) {
            // left coundary condition
            auto denominator = knots[i+degree]-knots[i];
            if (_floatIsZero(denominator)) {
                derCoeffs[i] = 0.0f;
            } else {
                derCoeffs[i] = degree*coeffs[0]/denominator;
            }
        } else if (i == n) {
            // right coundary condition
            auto denominator = knots[i+degree]-knots[i];
            if (_floatIsZero(denominator)) {
                derCoeffs[i] = 0.0f;
            } else {
                derCoeffs[i] = -degree*coeffs[i-1]/denominator;
            }
        } else {
            auto denominator = knots[i+degree]-knots[i];
            if (_floatIsZero(denominator)) {
                derCoeffs[i] = 0.0f;
            } else {
                derCoeffs[i] = degree*(coeffs[i]-coeffs[i-1])/denominator; 
            }
        }
    }
    return derCoeffs;
}

// Compute the updated B-spline coefficients for a refined knot vector
// obtained by adding one new knot to an existsing.
// This is done using Böhm's method.
// NOTE: The new knot z must be within the knot vector.
void bohmsMethod(int degree,
                 float z,
                 const std::vector<float>& knotsIn,
                 const std::vector<float>& coeffsIn,
                 std::vector<float>& knotsOut,
                 std::vector<float>& coeffsOut) {
    int numKnots = knotsIn.size();
    int numCoeffs = coeffsIn.size();
    if (z < knotsIn[0] || z > knotsIn[numKnots-1]) throw std::runtime_error("New knot not within knot vector");
    int mu = getMu(z, knotsIn);
    
    // Create the refined knot vector
    knotsOut.resize(numKnots+1);
    for (int i = 0; i <= mu; i++) {
        knotsOut[i] = knotsIn[i];
    }
    knotsOut[mu+1] = z;
    for (int i = mu+2; i < numKnots+1; i++) {
        knotsOut[i] = knotsIn[i-1];
    }

    // Compute the coefficients relative to refined knot vector
    coeffsOut.resize(numCoeffs+1);
    for (int i = 0; i <= mu-degree; i++) {
        coeffsOut[i] = coeffsIn[i];
    }
    for (int i = mu-degree+1; i <= mu; i++) {
        coeffsOut[i] = coeffsIn[i]*_specialDiv(z-knotsIn[i], knotsIn[i+degree]-knotsIn[i])
                        + coeffsIn[i-1]*_specialDiv(knotsIn[i+degree]-z, knotsIn[i+degree]-knotsIn[i]);
    }
    for (int i = mu+1; i < numCoeffs+1; i++) {
        coeffsOut[i] = coeffsIn[i-1];
    }

}

}   // end namespace bspline_storve
