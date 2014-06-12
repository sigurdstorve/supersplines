#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "bspline.hpp"

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