//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <vector>
#include <utility>  // for std::pair
#include <boost/test/unit_test.hpp>
#include "bspline.hpp"

// Workaround for MSVC2012 not supporting curly brace initialization for std::vector
template <typename T, size_t N>
std::vector<T> makeVector(const T(&data)[N]) {
    return std::vector<T>(data, data + N);
}

// Verfify that the simple implementation of linspace() works
// as expected. (Needed e.g. to generate not vectors)
BOOST_AUTO_TEST_CASE(LinspaceSanityCheck) {
    std::vector<float> v;
    bspline_storve::linspace(0.0f, 1.0f, 3, v);
    BOOST_CHECK_CLOSE(v[0], 0.0, 0.0001);
    BOOST_CHECK_CLOSE(v[1], 0.5, 0.0001);
    BOOST_CHECK_CLOSE(v[2], 1.0, 0.0001);

    bspline_storve::linspace(0.0f, 1.0f, 0, v);
    BOOST_CHECK_EQUAL(v.size(), 0);
    
    bspline_storve::linspace(0.0f, 10.0f, 1, v);
    BOOST_CHECK_CLOSE(v[0], 0.0, 0.0001);
    
    bspline_storve::linspace(0.0f, 10.0f, 2, v);
    BOOST_CHECK_CLOSE(v[0], 0.0, 0.0001);
    BOOST_CHECK_CLOSE(v[1], 10.0, 0.0001);
}

// Verify that the code tto create uniform regular knot vectors works.
BOOST_AUTO_TEST_CASE(UniformRegularKnotVectorTest) {
    int degree = 3;
    int numPoints = 7;
    float tStart = 1.0f;
    float tEnd = 2.0f;
    
    const float tempDesired[] = {1.0f, 1.0f, 1.0f, 1.0f, 1.25f, 1.50f, 1.75f, 2.0f, 2.0f, 2.0f, 2.0};
    std::vector<float> desired = makeVector(tempDesired);
    auto knots = bspline_storve::uniformRegularKnotVector(numPoints, degree, tStart, tEnd);
    BOOST_CHECK_EQUAL(knots.size(), 11);
    for (size_t i = 0; i < knots.size(); i++) {
        BOOST_CHECK_CLOSE(knots[i], desired[i], 0.0001);
    }
    
    // Make sure end-hack is not turned on by error
    int numKnots = knots.size();
    BOOST_CHECK(knots[numKnots-2] == knots[numKnots-1]);

    // Get the corresponding knot vector WITH end-hack
    knots = bspline_storve::uniformRegularKnotVector(numPoints, degree, tStart, tEnd, true);
    BOOST_CHECK(knots[numKnots-2] < knots[numKnots-1]);
    
}

// Verify that the analytical and recursive computations of B-splines
// give the same result when degree is 1.
BOOST_AUTO_TEST_CASE(CompareAnalyticalAndRecursiveDegree1) {
    int degree = 1;
    // Knots for 6 degree-1 B-splines.
    const float tempKnots[] = {-1.0f, 1.0f, 2.0f, 4.0f, 5.0f, 7.0f, 10.0f, 11.0f};
    std::vector<float> knots = makeVector(tempKnots);
    for (int basisFuncNo = 0; basisFuncNo < 6; basisFuncNo++) {
        for (float t = -2.0f; t <= 12.0f; t+=0.01f) {
            float recursive = bspline_storve::bsplineBasis(basisFuncNo,
                                                           degree,
                                                           t, knots);
            float analytic = bspline_storve::B1(basisFuncNo,
                                                t, knots);
            BOOST_CHECK_CLOSE(recursive, analytic, 0.0001);
        }
    }
}

// Verify that the analytical and recursive computations of B-splines
// give the same result when degree is 2.
BOOST_AUTO_TEST_CASE(CompareAnalyticalAndRecursiveDegree2) {
    int degree = 2;
    // Knots for 5 degree-2 B-splines.
    const float tempKnots[] = {-1.0f, 1.0f, 2.0f, 4.0f, 5.0f, 7.0f, 10.0f, 11.0f};
    std::vector<float> knots = makeVector(tempKnots);
    for (int basisFuncNo = 0; basisFuncNo < 5; basisFuncNo++) {
        for (float t = -2.0f; t <= 12.0f; t+=0.01f) {
            float recursive = bspline_storve::bsplineBasis(basisFuncNo,
                                                           degree,
                                                           t, knots);
            float analytic = bspline_storve::B2(basisFuncNo,
                                                t, knots);
            BOOST_CHECK_CLOSE(recursive, analytic, 0.0001);
        }
    }
}

// Verify that the analytical and recursive computations of B-splines
// give the same result when degree is 3.
BOOST_AUTO_TEST_CASE(CompareAnalyticalAndRecursiveDegree3) {
    int degree = 3;
    // Knots for 4 degree-3 B-splines.
    const float tempKnots[] = {-1.0f, 1.0f, 2.0f, 4.0f, 5.0f, 7.0f, 10.0f, 11.0f};
    std::vector<float> knots = makeVector(tempKnots);
    for (int basisFuncNo = 0; basisFuncNo < 4; basisFuncNo++) {
        for (float t = -2.0f; t <= 12.0f; t+=0.01f) {
            float recursive = bspline_storve::bsplineBasis(basisFuncNo,
                                                           degree,
                                                           t, knots);
            float analytic = bspline_storve::B3(basisFuncNo,
                                                t, knots);
            BOOST_CHECK_CLOSE(recursive, analytic, 0.0001);
        }
    }
}

// Verify that the knot span finder is correct
BOOST_AUTO_TEST_CASE(TestKnotSpanFinder) {
    const float tempKnots[] = {-1.0f, 1.0f, 2.0f, 4.0f, 5.0f, 7.0f, 10.0f, 11.0f};
    std::vector<float> knots = makeVector(tempKnots);
    std::vector<std::pair<float, int> > testPairs;
    testPairs.push_back(std::pair<float, int>(-1.0f, 0));
    testPairs.push_back(std::pair<float, int>(-0.5f, 0));
    testPairs.push_back(std::pair<float, int>(-0.01f, 0));
    testPairs.push_back(std::pair<float, int>(0.0f, 0));
    testPairs.push_back(std::pair<float, int>(0.5f, 0));
    testPairs.push_back(std::pair<float, int>(0.9f,  0));
    testPairs.push_back(std::pair<float, int>(1.0f, 1));
    testPairs.push_back(std::pair<float, int>(1.1f, 1));
    testPairs.push_back(std::pair<float, int>(1.99f, 1));
    testPairs.push_back(std::pair<float, int>(2.0f, 2));
    testPairs.push_back(std::pair<float, int>(3.4f, 2));
    testPairs.push_back(std::pair<float, int>(3.99f, 2));
    testPairs.push_back(std::pair<float, int>(4.0f, 3));
    testPairs.push_back(std::pair<float, int>(4.1f, 3));
    testPairs.push_back(std::pair<float, int>(4.99f, 3));
    testPairs.push_back(std::pair<float, int>(5.0f, 4));
    testPairs.push_back(std::pair<float, int>(6.0f, 4));
    testPairs.push_back(std::pair<float, int>(6.99f, 4));
    testPairs.push_back(std::pair<float, int>(7.0f, 5));
    testPairs.push_back(std::pair<float, int>(9.0f, 5));
    testPairs.push_back(std::pair<float, int>(9.99f, 5));
    testPairs.push_back(std::pair<float, int>(10.0f, 6));
    testPairs.push_back(std::pair<float, int>(10.5f, 6));
    testPairs.push_back(std::pair<float, int>(10.99f, 6));
        
    int lastMu = 1000; // crazy guess should not cause error
    for (auto it = testPairs.begin(); it != testPairs.end(); ++it) {
        float t = it->first;
        int desiredMu = it->second;
        
        // Without hints
        BOOST_CHECK_EQUAL(bspline_storve::getMu(t, knots), desiredMu);
        
        // With hints
        BOOST_CHECK_EQUAL(bspline_storve::getMu(t, knots, lastMu), desiredMu);
    }
}

// Verfify that the dimensionality of B-spline matrices are correct
BOOST_AUTO_TEST_CASE(TestBSplineMatrixDimensions) {
    int mu = 0;
    float x = 0.0f;
    for (int k = 1; k < 8; k++) {
        std::vector<float> knots(10*k);
        Eigen::MatrixXf m = bspline_storve::bsplineMatrix(k, mu, x, knots);
        // B-spline matrix k should be of dimensions k x (k+1)
        BOOST_CHECK_EQUAL(m.rows(), k);
        BOOST_CHECK_EQUAL(m.cols(), k+1);
    }

};

// Verify that the control polygon for a p+1-regular knot
// vector for n control points also has length n
BOOST_AUTO_TEST_CASE(TestControlPolygonAndRegularKnorVector) {
    float tStart = -1.0f;
    float tEnd = 2.2f;
    for (int degree = 0; degree < 10; degree++) {
        for (int n = 2*(degree+1); n < 100; n++) {
            // Regular            
            std::vector<float> knots;
            std::vector<float> controlPolygonTimes;
            knots = bspline_storve::uniformRegularKnotVector(n, degree, tStart, tEnd);
            controlPolygonTimes = bspline_storve::controlPolygon(degree, knots);
            BOOST_CHECK_EQUAL(controlPolygonTimes.size(), n);

            // With end hack
            knots = bspline_storve::uniformRegularKnotVector(n, degree, tStart, tEnd);
            controlPolygonTimes = bspline_storve::controlPolygon(degree, knots);
            BOOST_CHECK_EQUAL(controlPolygonTimes.size(), n);
        }
    }
}

// Verify correctness of control polygon times in some concrete cases
BOOST_AUTO_TEST_CASE(TestControlPolygonConcreteCases) {
    const float tempKnots[] = {1.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    std::vector<float> knots = makeVector(tempKnots);
    std::vector<float> cpTimes;
    int degree;

    degree = 1;
    cpTimes = bspline_storve::controlPolygon(degree, knots);
    BOOST_CHECK_CLOSE(cpTimes[0], knots[1], 0.0001);
    BOOST_CHECK_CLOSE(cpTimes[1], knots[2], 0.0001);
    BOOST_CHECK_CLOSE(cpTimes[2], knots[3], 0.0001);
    
    degree = 2;
    cpTimes = bspline_storve::controlPolygon(degree, knots);
    BOOST_CHECK_CLOSE(cpTimes[0], (knots[1]+knots[2])/degree, 0.0001);
    BOOST_CHECK_CLOSE(cpTimes[1], (knots[2]+knots[3])/degree, 0.0001);

    degree = 3;
    cpTimes = bspline_storve::controlPolygon(degree, knots);
    BOOST_CHECK_CLOSE(cpTimes[0], (knots[1]+knots[2]+knots[3])/degree, 0.0001);
}

// Verify that the cubic B-splines are correctly computed as products
// of B-spline matrices.
BOOST_AUTO_TEST_CASE(TestCubicBSplineMatrixProduct) {
    const float tempKnots[] = {-2.0f, -0.5f, 0.0f, 1.0f, 2.0f, 2.5f, 3.3f, 4.0f, 5.0f, 6.0f, 7.5f, 8.0f, 9.0f, 10.0f};
    std::vector<float> knots = makeVector(tempKnots);
    int degree = 3;
    int lastMu = 1000; // crazy guess should not cause errors
    for (float t = 1.0f; t < 5.99f; t += 0.01f) {
        int mu = bspline_storve::getMu(t, knots, lastMu);
        Eigen::VectorXf bsplineVec = bspline_storve::bsplineVector(degree, mu, t, knots);

        // Verify that we got degree+1 B-splines at t
        BOOST_CHECK_EQUAL(bsplineVec.size(), degree+1);
        
        // Verify that each element of B-spline vector is
        // equal to corresponding recursively computed value
        for (int i = 0; i < bsplineVec.size(); i++) {
            int basisFuncNo = mu-degree+i;
            float recursive = bspline_storve::bsplineBasis(basisFuncNo, degree, t, knots);
            BOOST_CHECK_CLOSE(bsplineVec[i], recursive, 0.0001);
        }
    }

}

// Verify that the computation of B-spline coefficients for the
// derivative of a given B-spline is correct for a simple degree-1 case
BOOST_AUTO_TEST_CASE(TestSplineDerivativeCoeffsComputationDegree1) {
    int degree = 1;
    std::vector<float> knots;
    std::vector<float> coeffs;
    coeffs.push_back(1.0f);
    
    // Test case: uniform knot spacing
    const float tempKnots1[] = {1.0f, 2.0f, 3.0f};
    knots = makeVector(tempKnots1);
    std::vector<float> derCoeffs = bspline_storve::computeDerivativeCoeffs(degree, coeffs, knots);
    BOOST_CHECK_EQUAL(derCoeffs.size(), 2);
    BOOST_CHECK_CLOSE(derCoeffs[0], 1.0, 0.0001);
    BOOST_CHECK_CLOSE(derCoeffs[1], -1.0, 0.0001);

    // Test case: non-uniform spacing
    const float tempKnots2[] = {-1.0f, 0.0f, 2.0f};
    knots = makeVector(tempKnots2);
    derCoeffs = bspline_storve::computeDerivativeCoeffs(degree, coeffs, knots);
    BOOST_CHECK_EQUAL(derCoeffs.size(), 2);
    BOOST_CHECK_CLOSE(derCoeffs[0], 1.0, 0.0001);
    BOOST_CHECK_CLOSE(derCoeffs[1], -0.5, 0.0001);

    // Test case: non-uniform spacing and more than one coeff
    const float tempCoeffs[] = {0.0f, 1.0f, 0.0f};
    coeffs = makeVector(tempCoeffs);
    const float tempCoeffs3[] = {1.0f, 2.0f, 4.0f, 8.0f, 16.0f};
    knots = makeVector(tempCoeffs3);
    derCoeffs = bspline_storve::computeDerivativeCoeffs(degree, coeffs, knots);
    BOOST_CHECK_EQUAL(derCoeffs.size(), 4);
    BOOST_CHECK_CLOSE(derCoeffs[0], 0.0, 0.0001);
    BOOST_CHECK_CLOSE(derCoeffs[1], 1.0/(knots[2]-knots[1]), 0.0001);
    BOOST_CHECK_CLOSE(derCoeffs[2], -1.0/(knots[3]-knots[2]), 0.0001);
    BOOST_CHECK_CLOSE(derCoeffs[3], 0.0, 0.0001);

}

// Simple test of B�hm's method
BOOST_AUTO_TEST_CASE(TestBohmsMethodSimple) {
    int degree;
    std::vector<float> knots, coeffs;
    std::vector<float> knotsOut, coeffsOut;

    degree = 0;
    const float tempKnots[] = {1.0f, 2.0f, 3.0f, 4.0f};
    knots = makeVector(tempKnots);
    const float tempCoeffs[] = {1.0f, -1.0, 2.0f};
    coeffs = makeVector(tempCoeffs);
    float newKnot = 2.5f;
    bspline_storve::bohmsMethod(degree, newKnot, knots, coeffs, knotsOut, coeffsOut);
    BOOST_REQUIRE_CLOSE(knotsOut[0], 1.0f, 0.0001);
    BOOST_REQUIRE_CLOSE(knotsOut[1], 2.0f, 0.0001);
    BOOST_REQUIRE_CLOSE(knotsOut[2], 2.5f, 0.0001);
    BOOST_REQUIRE_CLOSE(knotsOut[3], 3.0f, 0.0001);
    BOOST_REQUIRE_CLOSE(knotsOut[4], 4.0f, 0.0001);
    BOOST_REQUIRE_CLOSE(coeffsOut[0], 1.0f, 0.0001);
    BOOST_REQUIRE_CLOSE(coeffsOut[1], -1.0f, 0.0001);
    BOOST_REQUIRE_CLOSE(coeffsOut[2], -1.0f, 0.0001);
    BOOST_REQUIRE_CLOSE(coeffsOut[3], 2.0f, 0.0001);
}


// Test that a spline is unchanged when using B�hm's metod for knot insertion
BOOST_AUTO_TEST_CASE(TestBohmsMethodUnchagedSpline) {
    float tStart = 1.0f;
    float tEnd = 2.0f;
    int numCs = 25;
    int numCheckPoints = 50;

    // One big vector of coefficients - used by all tests
    std::vector<float> originalCoeffs;
    for (int i = 0; i < numCs; i++) {
        originalCoeffs.push_back(static_cast<float>(-10 + i*i));
    }

    // Generate the vector of parameter values where the original and
    // refined should be compared
    std::vector<float> checkParValues;
    bspline_storve::linspace(tStart-1, tEnd+1, numCheckPoints, checkParValues);

    for (int degree = 0; degree < 5; degree++) {
        std::vector<float> originalKnots = bspline_storve::uniformRegularKnotVector(numCs, degree, tStart, tEnd);
        std::vector<float> knots = originalKnots;
        std::vector<float> coeffs = originalCoeffs;

        // Add a knot and recompute coeffs
        for (float newKnot = tStart; newKnot < tEnd; newKnot += 0.01f) {
            std::vector<float> knotsOut, coeffsOut;
            bspline_storve::bohmsMethod(degree, newKnot, knots, coeffs, knotsOut, coeffsOut);
            BOOST_CHECK_EQUAL(knotsOut.size(), knots.size()+1);
            BOOST_CHECK_EQUAL(coeffsOut.size(), coeffs.size()+1);
            // Continue refining the refined.
            knots = knotsOut;
            coeffs = coeffsOut;

            // Evaluate and compare the original and cumulatively refined - they shall always be the same.
            auto originalEval = bspline_storve::renderSpline(degree, originalKnots, originalCoeffs, checkParValues);
            auto refinedEval = bspline_storve::renderSpline(degree, knotsOut, coeffsOut, checkParValues);
            for (size_t i = 0; i < originalEval.size(); i++) {
                BOOST_CHECK_CLOSE(originalEval[i], refinedEval[i], 0.0001);
            }
        }
    }

}
/*
BOOST_AUTO_TEST_CASE(TestZeroDetectionWhenNoZerosExists) {
    std::vector<float> knots = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f};
    std::vector<float> coeffs = {2.0f, 4.0f, 8.0f};
    int degree = 3;
    int numIterations = 10;
    auto zeroTimes = bspline_storve::computeZeros(degree, knots, coeffs, numIterations);
    BOOST_CHECK_EQUAL(zeroTimes.size(), 0);    
}
*/
BOOST_AUTO_TEST_CASE(TestZeroDetection1) {
    float tStart = 0.0f;
    float tEnd = 1.0f;
    int degree = 3;
    const float tempCoeffs[] = {0.1f, 1.0f, -8.0f, -0.1f};
    std::vector<float> coeffs = makeVector(tempCoeffs);
    std::vector<float> knots = bspline_storve::uniformRegularKnotVector(coeffs.size(), degree, tStart, tEnd);
    
    int numIterations = 10;
    auto zeroTimes = bspline_storve::computeZeros(degree, knots, coeffs, numIterations);
    BOOST_CHECK_EQUAL(zeroTimes.size(), 1);    
}