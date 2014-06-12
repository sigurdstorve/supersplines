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


int main(int argc, char** argv) {
    QApplication app(argc, argv);

    
    QWidget* mainWidget = new QWidget;
    mainWidget->resize(1024, 768);
    QHBoxLayout* hlayout = new QHBoxLayout;
    mainWidget->setLayout(hlayout);


    // Create a curve to plot
    std::vector<float> knots = {0.0, 0.0, 0.0, 0.0, 0.33, 0.67, 1.0, 1.0, 1.0, 1.0};
    std::vector<float> controlPoints = {0.0, 1.0, 2.0, -1.0, 0.0, 1.0};
    float x = 0.5;
    int p = 3;
    std::vector<float> xs;
    std::vector<float> ys;
    for (float t = 0; t < 1.0; t+=0.001) {
        xs.push_back(t);
        //ys.push_back(bspline_storve::bsplineBasis(0, p, t, knots));
    }
    ys = bspline_storve::renderSpline(p, knots, controlPoints, xs);

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
    std::vector<float> cpTimes = bspline_storve::controlPolygon(p, knots);
    auto cpTimesDouble = stdVectorToDouble(cpTimes);
    auto cpValuesDouble = stdVectorToDouble(controlPoints);
    curve2->setSamples(cpTimesDouble.data(), cpValuesDouble.data(), cpTimesDouble.size());
    curve2->attach(plot);
    
    mainWidget->show();
    std::cout << "Qt version: ";
    std::cout << QT_VERSION_STR << std::endl;
    
    return app.exec();
}
