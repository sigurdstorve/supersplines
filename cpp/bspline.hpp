#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <stdexcept>

// TODO: Stuff using Eigen3 is hardcoded to use float and not
// double. Is this ok?

namespace bspline_storve {

// Simple implementation of functionality like NumPy's linspace()
template <typename T>
void linspace(T left, T right, int n, std::vector<T>& res) {
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
template <typename T>
T bsplineBasis(int j, int p, T x, const std::vector<T>& knots) {
    if (p == 0 ) {
        if ( (x >= knots[j]) && (x < knots[j+1]) ) {
            return 1.0;
        } else {
            return 0.0;
        }
    } else {
        T left = _specialDiv((x-knots[j])*bsplineBasis(j,p-1,x,knots), knots[j+p]-knots[j]);
        T right = _specialDiv((knots[j+1+p]-x)*bsplineBasis(j+1,p-1,x,knots), knots[j+1+p]-knots[j+1]);
        return left + right;
    }
}

/*
    Find the knot interval index (mu) for a given value x
    This is a prerequisite for alg 2.20 and 2.21.
    Optional to supply a suggestion to try first.
    Throws on error.
*/
template <typename T>
int getMu(T x, const std::vector<T>& knots, int mu=-1) {
    if (mu != -1 && (mu <= static_cast<int>(knots.size())-2) && (mu >= 0)) {
        if (x >= knots[mu] && x < knots[mu+1]) {
            return mu;
        }
    }
    for (int mu = 0; mu < (static_cast<int>(knots.size())-1); mu++) {
        if (x >= knots[mu] && x < knots[mu+1]) {
            return mu;
        }
    }
    throw std::runtime_error("Illegal knot vector or x value");
}

// Compute B-spline matrix from theorem 2.18
template <typename T>
Eigen::MatrixXf bsplineMatrix(int k, int mu, T x, const std::vector<T>& knots) {
    if (k <= 0) throw std::runtime_error("k must be greater than zero");
    Eigen::MatrixXf res(k, k+1);
    for (int row = 0; row <k; row++) {
        for (int col = 0; col < k+1; col++) {
            res(row,col) = 0.0f;
        }
    }    
    for (int row = 0; row < k; row++) {
        T common_denom = knots[mu+(row+1)] - knots[mu+(row+1)-k];
        res(row,row) = _specialDiv(knots[mu+(row+1)] - x, common_denom);
        res(row,row+1) = _specialDiv(x - knots[mu + (row+1) - k], common_denom);
    }
    return res;
}

// Compute B-spline vector using algorithm 2.19/2.20
// Returns (B_{\mu-p,p}(x),...,B_{\mu,p}(x))
template <typename T>
Eigen::VectorXf bsplineVector(int p, int mu, T x, const std::vector<T>& knots) {
    Eigen::MatrixXf res = bsplineMatrix(1, mu, x, knots);
    for (int k = 2; k <= p; k++) {
        res *= bsplineMatrix(k, mu, x, knots);
    }
    return res;
}

// Direkte
template <typename T>
T B0(int j, T x, const std::vector<T>& knots) {
    if (knots[j] <= x && x < knots[j+1]) {
        return 1.0;
    } else {
        return 0.0;
    }
}

template <typename T>
T B1(int j, T x, const std::vector<T>& knots) {
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
        return 0.0;
    } 
}

template <typename T>
T B2(int j, T x, const std::vector<T>& knots) {
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
        return 0.0;
    }
}

template <typename T>
T B3(int j, T x, const std::vector<T>& knots) {
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
        return 0.0;
    }
}

/*
Create a p+1-regular uniform knot vector for
a given number of control points
Throws if n is too small
endHack: Make the last knot bigger than the one before it, which has the
effect of making it non-zero in last knot value.
*/
template <typename T>
std::vector<T> uniformRegularKnotVector(int numPoints,
                                        int degree,
                                        T tStart=0.0,
                                        T tEnd=1.0,
                                        bool endHack=false) {

    // The minimum length of a degree+1-regular knot vector is 2*(degree+1)
    if (numPoints < degree+1) {
        throw std::runtime_error("Too small n for a uniform regular knot vector");
    }
    
    std::vector<T> knots;
    // degree+1 copies of tStart left and degree+1 copies of tEnd right
    // but one of each in linspace
    for (int i = 0; i < degree; i++) {
        knots.push_back(tStart);
    }
    std::vector<T> temp;
    linspace(tStart, tEnd, numPoints+1-degree, temp);
    for (T t : temp) {
        knots.push_back(t);
    }
    for (int i = 0; i < degree; i++) {
        knots.push_back(tEnd);
    }
    
    if (endHack) {
        int numKnots = knots.size();
        knots[numKnots-1] += 1.0;
    }
    return knots;
}

/*
Return the control point abscissa for the control polygon
of a one-dimensional spline.
*/
template <typename T>
std::vector<T> controlPolygon(int p, const std::vector<T>& knots) {
    int n = knots.size() - p - 1;
    std::vector<T> abscissas;
    for (int i = 0; i < n; i++) { 
        T sum = 0.0;
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
template <typename T>
std::vector<T> renderSpline(int p,
                            const std::vector<T>& knots,
                            const std::vector<T>& controlPoints,
                            const std::vector<T>& ts) {
    std::vector<T> ys;
    for (T t : ts) {
        T y = 0.0; 
        for (size_t j = 0; j < controlPoints.size(); j++) {
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
template <typename T>
std::vector<T> leastSquaresFit(const std::vector<T>& xs,
                               const std::vector<T>& ys,
                               const std::vector<T>& knots,
                               int p,
                               const std::vector<T>& _weights = std::vector<T>()) {
    int numKnots = knots.size();
    int numDataPoints = xs.size();
    //assert(ys.size() == numDataPoints);
   
    // The number of knots determines the number of control points
    // to use in the approximation.
    int numApproxPoints = numKnots - p - 1;
   
    std::vector<T> weights(numDataPoints);
    if (_weights.size() == 0) {
        for (int i = 0; i < numDataPoints; i++) weights[i] = 1.0;
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
    std::vector<T> res(numApproxPoints);
    for (int i = 0; i < numApproxPoints; i++) {
        res[i] = temp(i);
    }
    return res;
}

// Compute the control points of the derivative of a given spline.
// The derivatoive of a degree-p spline with n B-spline coefficients
// is a degree-(p-1) spline with n+1 control points and can be
// expressed on the same knot vector.
template <typename T>
std::vector<T> computeDerivativeCoeffs(int degree,
                                       const std::vector<T>& coeffs,
                                       const std::vector<T>& knots) {
    if (degree <= 0) throw std::runtime_error("Cannot differentiate a spline of degree lower than 1");
    int n = coeffs.size();
    std::vector<T> derCoeffs(n+1);
    // Handle boundary cases
    for (int i = 0; i < n+1; i++) {
        if (i == 0) {
            // left coundary condition
            auto denominator = knots[i+degree]-knots[i];
            if (_floatIsZero(denominator)) {
                derCoeffs[i] = 0.0;
            } else {
                derCoeffs[i] = degree*coeffs[0]/denominator;
            }
        } else if (i == n) {
            // right coundary condition
            auto denominator = knots[i+degree]-knots[i];
            if (_floatIsZero(denominator)) {
                derCoeffs[i] = 0.0;
            } else {
                derCoeffs[i] = -degree*coeffs[i-1]/denominator;
            }
        } else {
            auto denominator = knots[i+degree]-knots[i];
            if (_floatIsZero(denominator)) {
                derCoeffs[i] = 0.0;
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
template <typename T>
void bohmsMethod(int degree,
                 T z,
                 const std::vector<T>& knotsIn,
                 const std::vector<T>& coeffsIn,
                 std::vector<T>& knotsOut,
                 std::vector<T>& coeffsOut) {
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

// Experimental. Does not take identically equal to zero on an interval
// possibility into account.
// TODO: Fix handling of end-points and if a zero-crossing IS the zero exactly
// (currently this breaks the algorithm)
template <typename T>
std::vector<T> computeZeros(int degree,
                            std::vector<T> /*not reference!*/ knots,
                            std::vector<T> /*not reference!*/ coeffs,
                            int numIterations) {
    std::vector<T> zeroTimes;
    for (int iteration = 0; iteration < numIterations; iteration++) {
        int numCoeffs = coeffs.size();
        // Get the (updated) control polygon times
        auto cpTimes = controlPolygon(degree, knots);
        std::vector<T> tempTimes;
        //std::cout << "*** Start of a new iteration\n";
        //std::cout << "Coeffs:\n";
        //for (int i = 0; i < numCoeffs; i++) {
        //    std::cout << coeffs[i] << std::endl;
        //}
        for (int i = 1; i < numCoeffs; i++) {
            //std::cout << coeffs[i-1]*coeffs[i] << std::endl;
            if (coeffs[i-1]*coeffs[i] < 0.0) {
                // Compute exactly where the zero-crossings occurs.
                T t = (cpTimes[i-1]*coeffs[i]-cpTimes[i]*coeffs[i-1])/(coeffs[i]-coeffs[i-1]);
                //std::cout << "zero at t=" << t << std::endl;
                tempTimes.push_back(t);
            }
        }
        // Insert a new knot at every zero-crossing that was found.
        for (T zeroCrossTime : tempTimes) {            
            std::vector<T> knotsOut, coeffsOut;
            //std::cout << "inserting knot at " << zeroCrossTime << std::endl;
            bohmsMethod(degree, zeroCrossTime, knots, coeffs, knotsOut, coeffsOut);
            knots = knotsOut;
            coeffs = coeffsOut;
        }
        //std::cout << "Found " << tempTimes.size() << " zero crossings\n";
        //std::cout << "Current knot vector length: " << knots.size() << std::endl;
        zeroTimes = tempTimes;
    }
    return zeroTimes;
}

}   // end namespace bspline_storve
