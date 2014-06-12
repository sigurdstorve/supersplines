#include <QApplication>
#include <QPen>
#include <QWidget>
#include <QHBoxLayout>
#ifdef _WIN32
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#else
#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#endif
#include <vector>
#include <cstdlib>
#include "bspline.hpp"

// A simple Qt application containing a single QwtPlot to plot a curve

std::vector<float> stdVectorToFloat(const std::vector<double>& in) {
    std::vector<float> res(in.begin(), in.end());
    return res;
}

std::vector<double> stdVectorToDouble(const std::vector<float>& in) {
    std::vector<double> res(in.begin(), in.end());
    return res;
}


// Show control polygon and rendered spline function.
int main1(int argc, char** argv) {
    QApplication app(argc, argv);

    
    QWidget* mainWidget = new QWidget;
    mainWidget->resize(1024, 768);
    QHBoxLayout* hlayout = new QHBoxLayout;
    mainWidget->setLayout(hlayout);

    int degree = 3; // cubic

    // Create a curve to plot
    std::vector<float> controlPoints = {0.0, 1.0, 2.0, -1.0, 0.0, 1.0};
    std::vector<float> knots = bspline_storve::uniformRegularKnotVector(controlPoints.size(),
                                                                        degree,
                                                                        0.0f, 1.0f, true);
    std::vector<float> xs;
    std::vector<float> ys;
    for (float t = 0; t < 1.0; t+=0.001) {
        xs.push_back(t);
    }
    ys = bspline_storve::renderSpline(degree, knots, controlPoints, xs);

    QwtPlot* plot = new QwtPlot;
    hlayout->addWidget(plot);
    QwtPlotCurve* curve1 = new QwtPlotCurve("Curve 1");
    curve1->setPen(QPen(QColor(255, 0, 0)));
    auto xsDouble = stdVectorToDouble(xs);
    auto ysDouble = stdVectorToDouble(ys);
    curve1->setSamples(xsDouble.data(), ysDouble.data(), xsDouble.size());
    curve1->attach(plot);
    
    QwtPlotCurve* curve2 = new QwtPlotCurve("Control Polygon");
    curve2->setPen(QPen(QColor(0, 0, 0)));
    std::vector<float> cpTimes = bspline_storve::controlPolygon(degree, knots);
    auto cpTimesDouble = stdVectorToDouble(cpTimes);
    auto cpValuesDouble = stdVectorToDouble(controlPoints);
    curve2->setSamples(cpTimesDouble.data(), cpValuesDouble.data(), cpTimesDouble.size());
    curve2->attach(plot);
    
    mainWidget->show();
    std::cout << "Qt version: ";
    std::cout << QT_VERSION_STR << std::endl;
    
    return app.exec();
}

// Simple least squares spline function approximation demo.
int main(int argc, char** argv) {
    QApplication app(argc, argv);
    
    QWidget* mainWidget = new QWidget;
    mainWidget->resize(1024, 768);
    QHBoxLayout* hlayout = new QHBoxLayout;
    mainWidget->setLayout(hlayout);

    int degree = 3;
    int numSamples = 50;
    int numControlPoints = 12;
    
    // Generate some random points
    std::vector<float> dataXs;
    std::vector<float> dataYs;
    bspline_storve::linspace(0.0, 1.0, numSamples, dataXs);
    for (int i = 0; i < numSamples; i++) {
        dataYs.push_back(std::rand() % 100);
    }
    
    // Render the spline approximation.
    std::vector<float> knots = bspline_storve::uniformRegularKnotVector(numControlPoints,
                                                                        degree,
                                                                        dataXs[0],
                                                                        dataXs[numSamples-1],
                                                                        true);
                                                                        
    // Compute least squares spline approximation in spline space defined
    // by degree and knot vector.
    auto controlPoints = bspline_storve::leastSquaresFit(dataXs, dataYs, knots, degree);

    std::cout << "BAJS\n";

    std::vector<float> xs;
    std::vector<float> ys;
    bspline_storve::linspace(0.0, 1.0, 1000, xs);
    ys = bspline_storve::renderSpline(degree, knots, controlPoints, xs);

    QwtPlot* plot = new QwtPlot;
    hlayout->addWidget(plot);
    QwtPlotCurve* curve1 = new QwtPlotCurve("Curve 1");
    curve1->setPen(QPen(QColor(255, 0, 0)));
    auto xsDouble = stdVectorToDouble(xs);
    auto ysDouble = stdVectorToDouble(ys);
    curve1->setSamples(xsDouble.data(), ysDouble.data(), xsDouble.size());
    curve1->attach(plot);
    
    QwtPlotCurve* curve2 = new QwtPlotCurve("Control Polygon");
    curve2->setPen(QPen(QColor(0, 0, 0)));
    std::vector<float> cpTimes = bspline_storve::controlPolygon(degree, knots);
    auto cpTimesDouble = stdVectorToDouble(cpTimes);
    auto cpValuesDouble = stdVectorToDouble(controlPoints);
    curve2->setSamples(cpTimesDouble.data(), cpValuesDouble.data(), cpTimesDouble.size());
    curve2->attach(plot);
    
    QwtPlotCurve* curve3 = new QwtPlotCurve("Raw data");
    curve3->setPen(QPen(QColor(0,255,0)));
    auto dataXsDouble = stdVectorToDouble(dataXs);
    auto dataYsDouble = stdVectorToDouble(dataYs);
    curve3->setSamples(dataXsDouble.data(), dataYsDouble.data(), dataXsDouble.size());
    curve3->attach(plot);
    
    mainWidget->show();
    std::cout << "Qt version: ";
    std::cout << QT_VERSION_STR << std::endl;
    
    return app.exec();
}
