#include <vector>
#include <Eigen/Dense>
#include "bspline.hpp"

using namespace bspline_storve;

// Workaround for MSVC2012 not supporting curly brace initialization for std::vector
template <typename T, size_t N>
std::vector<T> makeVector(const T(&data)[N]) {
    return std::vector<T>(data, data + N);
}

std::vector<float> stdVectorToFloat(const std::vector<double>& in) {
    std::vector<float> res(in.begin(), in.end());
    return res;
}

void testBasisFunctions1() {
    // Enough knots for one cubic b-spline

    const float tempKnots[] = {0.0, 1.0, 1.0, 3.0, 4.0};
    std::vector<float> knots = makeVector(tempKnots);
    float x = 1.0;
    int p = 3;
    float b = bsplineBasis(0, p, x, knots);
    std::cout << "b = " << b << std::endl;
}

void testBSplineMatrix() {
    int k = 3;
    int mu = 3;
    const float tempKnots[] = {1,2,3,4,5,6,7,8,9,10,11,12};
    std::vector<float> knots = makeVector(tempKnots);
    float x = 1.0f;
    auto res = bsplineMatrix(k, mu, x, knots);
    std::cout << "Bk = " << res << std::endl;
}

void testDegree1() {
    std::cout << __FUNCTION__ << std::endl;
    const float tempKnots[] = {0.0, 1.0, 3.0};
    std::vector<float> knots = makeVector(tempKnots);
    bool pass = true;
    for (float x = -1.0; x < 4.0; x += 0.001f) {
        float b1 = bsplineBasis(0, 1, x, knots);
        float b2 = B1(0, x, knots);
        if (std::abs(b1-b2) > 1e-6) {
            std::cout << "Error: " << b1 << " " << b2 << std::endl;
            pass = false;
        }
    }
    if (!pass) {
        std::cout << __FUNCTION__ << " error\n";
    }
}

void testDegree2() {
    std::cout << __FUNCTION__ << std::endl;
    const float tempKnots[] = {0.0f, 1.0f, 3.0f, 4.0f};
    std::vector<float> knots = makeVector(tempKnots);
    bool pass = true;
    for (float x = -1.0; x < 5.0; x += 0.001f) {
        float b1 = bsplineBasis(0, 2, x, knots);
        float b2 = B2(0, x, knots);
        if (std::abs(b1-b2) > 1e-6) {
            std::cout << "Error: " << b1 << " " << b2 << std::endl;
            pass = false;
        }
    }
    if (!pass) {
        std::cout << __FUNCTION__ << " error\n";
    }
}

void testDegree3() {
    std::cout << __FUNCTION__ << std::endl;
    const float tempKnots[] = {0.0, 1.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<float> knots = makeVector(tempKnots);
    bool pass = true;
    for (int j = 0; j < 2; j++) {
        for (float x = -1.0; x < 7.0; x += 0.001f) {
            float b1 = bsplineBasis(j, 3, x, knots);
            float b2 = B3(j, x, knots);
            if (std::abs(b1-b2) > 1e-6) {
                std::cout << "Error: " << b1 << " " << b2 << std::endl;
                pass = false;
            }
        }
    }
    if (!pass) {
        std::cout << __FUNCTION__ << " error\n";
    }
}

void testLinAlg1() {
    std::cout << __FUNCTION__ << std::endl;
    Eigen::Matrix2d a;
    Eigen::MatrixXd b(2,2), c(2,2);
    a << 1,2,3,4;
    b << 2,3,1,4;
    a *= b;
   // std::cout << c << std::endl;
}

void testLinAlg2() {
 using namespace Eigen;
 Matrix2d a;
a << 1, 2,
3, 4;
Vector3d v(1,2,3);
std::cout << "a * 2.5 =\n" << a * 2.5 << std::endl;
std::cout << "0.1 * v =\n" << 0.1 * v << std::endl;
std::cout << "Doing v *= 2;" << std::endl;
v *= 2;
std::cout << "Now v =\n" << v << std::endl;
}

void testLeastSquares() {
    using namespace Eigen;
    MatrixXf A = MatrixXf::Random(3,2);
    std::cout << "A:" << A << std::endl;
    VectorXf b = VectorXf::Random(3);
    std::cout << "b: " << b << std::endl;
    std::cout << "LSQ solution: " << A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b) << std::endl;
}

void testLinAlg3() {
    using namespace Eigen;
    Matrix2d mat;
    mat << 1, 2,
    3, 4;
    Vector2d u(-1,1), v(2,0);
    std::cout << "Here is mat*mat:\n" << mat*mat << std::endl;
    std::cout << "Here is mat*u:\n" << mat*u << std::endl;
    std::cout << "Here is u^T*mat:\n" << u.transpose()*mat << std::endl;
    std::cout << "Here is u^T*v:\n" << u.transpose()*v << std::endl;
    std::cout << "Here is u*v^T:\n" << u*v.transpose() << std::endl;
    std::cout << "Let's multiply mat by itself" << std::endl;
    mat = mat*mat;
    std::cout << "Now mat is mat:\n" << mat << std::endl;
}

int main(void) {
    std::cout << "Test" << std::endl;
    
    testBasisFunctions1();
    testDegree1();
    testDegree2();
    testDegree3();
    testLinAlg3();
    testLeastSquares();
    testBSplineMatrix();

}
