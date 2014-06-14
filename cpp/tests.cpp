#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <vector>
#include <utility>  // for std::pair
#include <boost/test/unit_test.hpp>
#include "bspline.hpp"

// Verfify that the simple implementation of linspace() works
// as expected. (Needed e.g. to generate not vectors)
BOOST_AUTO_TEST_CASE(LinspaceSanityCheck) {
    std::vector<float> v;
    bspline_storve::linspace(0.0, 1.0, 3, v);
    BOOST_CHECK_CLOSE(v[0], 0.0, 0.0001);
    BOOST_CHECK_CLOSE(v[1], 0.5, 0.0001);
    BOOST_CHECK_CLOSE(v[2], 1.0, 0.0001);

    bspline_storve::linspace(0.0, 1.0, 0, v);
    BOOST_CHECK_EQUAL(v.size(), 0);
    
    bspline_storve::linspace(0.0, 10.0, 1, v);
    BOOST_CHECK_CLOSE(v[0], 0.0, 0.0001);
    
    bspline_storve::linspace(0.0, 10.0, 2, v);
    BOOST_CHECK_CLOSE(v[0], 0.0, 0.0001);
    BOOST_CHECK_CLOSE(v[1], 10.0, 0.0001);
}

// Verify that the code tto create uniform regular knot vectors works.
BOOST_AUTO_TEST_CASE(UniformRegularKnotVectorTest) {
    int degree = 3;
    int numPoints = 7;
    float tStart = 1.0f;
    float tEnd = 2.0f;
    
    std::vector<float> desired = {1.0f, 1.0f, 1.0f, 1.0f, 1.25f, 1.50f, 1.75f, 2.0f, 2.0f, 2.0f, 2.0};
    auto knots = bspline_storve::uniformRegularKnotVector(numPoints, degree, tStart, tEnd);
    BOOST_CHECK_EQUAL(knots.size(), 11);
    for (int i = 0; i < knots.size(); i++) {
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
    std::vector<float> knots = {-1.0f, 1.0f, 2.0f, 4.0f, 5.0f, 7.0f, 10.0f, 11.0f};
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
    std::vector<float> knots = {-1.0f, 1.0f, 2.0f, 4.0f, 5.0f, 7.0f, 10.0f, 11.0f};
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
    std::vector<float> knots = {-1.0f, 1.0f, 2.0f, 4.0f, 5.0f, 7.0f, 10.0f, 11.0f};
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
    std::vector<float> knots = {-1.0f, 1.0f, 2.0f, 4.0f, 5.0f, 7.0f, 10.0f, 11.0f};
    std::vector<std::pair<float, int> > testPairs{
        {-1.0f, 0}, {-0.5f, 0}, {-0.01f, 0},
        {0.0f, 0}, {0.5f, 0}, {0.9f,  0},
        {1.0f, 1}, {1.1f, 1}, {1.99f, 1},
        {2.0f, 2}, {3.4f, 2}, {3.99f, 2},
        {4.0f, 3}, {4.1f, 3}, {4.99f, 3},
        {5.0f, 4}, {6.0f, 4}, {6.99f, 4},
        {7.0f, 5}, {9.0f, 5}, {9.99f, 5},
        {10.0f, 6}, {10.5f, 6}, {10.99f, 6}
    };
    
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
    int x = 0.0f;
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
    std::vector<float> knots = {1.0f, 3.0f, 4.0f, 5.0f, 6.0f};
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
    std::vector<float> knots = {-2.0f, -0.5f, 0.0f, 1.0f, 2.0f, 2.5f, 3.3f, 4.0f, 5.0f, 6.0f, 7.5f, 8.0f, 9.0f, 10.0f};
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
    std::vector<float> coeffs = {1.0f};
    
    // Test case: uniform knot spacing
    knots = {1.0f, 2.0f, 3.0f};
    std::vector<float> derCoeffs = bspline_storve::computeDerivativeCoeffs(degree, coeffs, knots);
    BOOST_CHECK_EQUAL(derCoeffs.size(), 2);
    BOOST_CHECK_CLOSE(derCoeffs[0], 1.0, 0.0001);
    BOOST_CHECK_CLOSE(derCoeffs[1], -1.0, 0.0001);

    // Test case: non-uniform spacing
    knots = {-1.0f, 0.0f, 2.0f};
    derCoeffs = bspline_storve::computeDerivativeCoeffs(degree, coeffs, knots);
    BOOST_CHECK_EQUAL(derCoeffs.size(), 2);
    BOOST_CHECK_CLOSE(derCoeffs[0], 1.0, 0.0001);
    BOOST_CHECK_CLOSE(derCoeffs[1], -0.5, 0.0001);

    // Test case: non-uniform spacing and more than one coeff
    coeffs = {0.0f, 1.0f, 0.0f};
    knots = {1.0f, 2.0f, 4.0f, 8.0f, 16.0f};
    derCoeffs = bspline_storve::computeDerivativeCoeffs(degree, coeffs, knots);
    BOOST_CHECK_EQUAL(derCoeffs.size(), 4);
    BOOST_CHECK_CLOSE(derCoeffs[0], 0.0, 0.0001);
    BOOST_CHECK_CLOSE(derCoeffs[1], 1.0/(knots[2]-knots[1]), 0.0001);
    BOOST_CHECK_CLOSE(derCoeffs[2], -1.0/(knots[3]-knots[2]), 0.0001);
    BOOST_CHECK_CLOSE(derCoeffs[3], 0.0, 0.0001);

}