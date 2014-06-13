#include <QPen>
#include <QWidget>
#include <QHBoxLayout>
#include <QSlider>
#include <QPushButton>
#ifdef _WIN32
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_symbol.h>
#else
#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>
#endif
#include <vector>
#include <algorithm>
#include "bspline.hpp"
#include "CsvHandler.hpp"
#include "CsvHandler.cpp" // HACK-O-RAMA

#include "MyMainWindow.hpp"

std::vector<float> stdVectorToFloat(const std::vector<double>& in) {
    std::vector<float> res(in.begin(), in.end());
    return res;
}

std::vector<double> stdVectorToDouble(const std::vector<float>& in) {
    std::vector<double> res(in.begin(), in.end());
    return res;
}

MyMainWindow::MyMainWindow(const std::string& csvFile, QWidget* parent) 
        : QMainWindow(parent) {
    QWidget* mainWidget = new QWidget;
    mainWidget->resize(1024, 768);
    QHBoxLayout* hlayout = new QHBoxLayout;
    mainWidget->setLayout(hlayout);

    m_slider = new QSlider;
    m_slider->setValue(m_numControlPoints);
    connect(m_slider, SIGNAL(valueChanged(int)),
            this, SLOT(setNumControlPoints(int)));
    hlayout->addWidget(m_slider);
    
    m_button = new QPushButton("Control polygon on/off");
    hlayout->addWidget(m_button);
    connect(m_button, SIGNAL(released()), this, SLOT(toggleControlPolygon()));
    
    m_degree = 3;
    m_numControlPoints = 20;
    m_showControlPolygon = true;

    // Prepare plotting
    m_plot = new QwtPlot;
    hlayout->addWidget(m_plot);
    m_curve1 = new QwtPlotCurve("Least squares spline approx");
    m_curve2 = new QwtPlotCurve("Control Polygon");
    m_curve3 = new QwtPlotCurve("Raw data");
    m_curve1->setRenderHint(QwtPlotItem::RenderAntialiased);
    auto sym = new QwtSymbol(QwtSymbol::XCross);
    m_curve1->attach(m_plot);
    m_curve2->setRenderHint(QwtPlotItem::RenderAntialiased);
    sym->setSize(QSize(20,20));
    m_curve2->setSymbol(sym);
    m_curve2->attach(m_plot);
    m_curve3->setRenderHint(QwtPlotItem::RenderAntialiased);
    m_curve3->attach(m_plot);
    
    loadCsvData(csvFile);
    updateSplineApprox();
    mainWidget->show();
}
void MyMainWindow::loadCsvData(const std::string& csvFile) {
    loadTimeSignal(csvFile, m_dataXs, m_dataYs);
    int numSamples = m_dataXs.size();

    // Update slider limits to reflect new data
    m_slider->setMinimum(m_degree+1);
    m_slider->setMaximum(numSamples-1);
    
    // Update xlim and ylim
    m_plot->setAxisScale(QwtPlot::xBottom, m_dataXs.front(), m_dataXs.back());
    
    float yMin = *std::min_element(m_dataYs.begin(), m_dataYs.end());
    float yMax = *std::max_element(m_dataYs.begin(), m_dataYs.end());
    float yLen = yMax-yMin;
    m_plot->setAxisScale(QwtPlot::yLeft, yMin-0.1*yLen, yMax+0.1*yLen);
}

void MyMainWindow::updateSplineApprox() {
    
    // Create uniform regular knot vector
    int numSamples = m_dataXs.size();
    std::vector<float> knots = bspline_storve::uniformRegularKnotVector(m_numControlPoints,
                                                                        m_degree,
                                                                        m_dataXs[0],
                                                                        m_dataXs[numSamples-1],
                                                                        true);
                                                                        
    // Compute least squares spline approximation in spline space defined
    // by degree and knot vector.
    auto controlPoints = bspline_storve::leastSquaresFit(m_dataXs, m_dataYs, knots, m_degree);

    // Render the spline approximation.
    std::vector<float> xs;
    std::vector<float> ys;
    bspline_storve::linspace(m_dataXs[0], m_dataXs[numSamples-1], 1000, xs);
    ys = bspline_storve::renderSpline(m_degree, knots, controlPoints, xs);

    m_curve1->setPen(QPen(QColor(255, 0, 0)));
    auto xsDouble = stdVectorToDouble(xs);
    auto ysDouble = stdVectorToDouble(ys);
    m_curve1->setSamples(xsDouble.data(), ysDouble.data(), xsDouble.size());
    
    m_curve2->setPen(QPen(QColor(0, 0, 0)));
    std::vector<float> cpTimes = bspline_storve::controlPolygon(m_degree, knots);
    auto cpTimesDouble = stdVectorToDouble(cpTimes);
    auto cpValuesDouble = stdVectorToDouble(controlPoints);
    m_curve2->setSamples(cpTimesDouble.data(), cpValuesDouble.data(), cpTimesDouble.size());
    
    m_curve3->setPen(QPen(QColor(0,0,127)));
    auto dataXsDouble = stdVectorToDouble(m_dataXs);
    auto dataYsDouble = stdVectorToDouble(m_dataYs);
    m_curve3->setSamples(dataXsDouble.data(), dataYsDouble.data(), dataXsDouble.size());
    
    m_curve2->setVisible(m_showControlPolygon);
    
    m_plot->replot();
}