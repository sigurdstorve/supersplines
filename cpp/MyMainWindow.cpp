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
#include "QTweakWindow.hpp"

#include "MyMainWindow.hpp"

std::vector<float> stdVectorToFloat(const std::vector<double>& in) {
    std::vector<float> res(in.begin(), in.end());
    return res;
}

std::vector<double> stdVectorToDouble(const std::vector<float>& in) {
    std::vector<double> res(in.begin(), in.end());
    return res;
}

void MyMainWindow::refresh() {
    m_tweakWin->pushAll();
    updateSplineApprox();
}


MyMainWindow::MyMainWindow(const std::string& csvFile, QWidget* parent) 
        : QMainWindow(parent) {
    QWidget* mainWidget = new QWidget;
    mainWidget->resize(1024, 768);
    QHBoxLayout* hlayout = new QHBoxLayout;
    mainWidget->setLayout(hlayout);

    m_degree = 2;
    m_numControlPoints = 20;
    m_showControlPolygon = true;
    m_showDerivative = false;
    m_derivativeScale = 0.01;
    m_showBaseline = true;
    m_showDerivPolygon = false;
    m_showSplineFit = true;
    m_showDerivZeros = false;

    m_tweakWin = new QTweakWindow;
    m_tweakWin->registerVariable("m_degree", &m_degree, 0, 10);
    m_tweakWin->registerVariable("m_numControlPoints", &m_numControlPoints, 10, 100);
    m_tweakWin->registerVariable("m_showControlPolygon", &m_showControlPolygon);
    m_tweakWin->registerVariable("m_showDerivative", &m_showDerivative);
    m_tweakWin->registerVariable("m_derivativeScale", &m_derivativeScale, 0.0f, 2.0f);
    m_tweakWin->registerVariable("m_showBaseline", &m_showBaseline);
    m_tweakWin->registerVariable("m_showDerivPolygon", &m_showDerivPolygon);
    m_tweakWin->registerVariable("m_showSplineFit", &m_showSplineFit);
    m_tweakWin->registerVariable("m_showDerivZeros", &m_showDerivZeros);
    connect(m_tweakWin, SIGNAL(changesAvailable()), this, SLOT(refresh()));
    hlayout->addWidget(m_tweakWin);
    
    // Prepare plotting
    m_plot = new QwtPlot;
    hlayout->addWidget(m_plot);
    m_curve1 = new QwtPlotCurve("Least squares spline approx");
    m_curve2 = new QwtPlotCurve("Control Polygon");
    m_curve3 = new QwtPlotCurve("Raw data");
    m_curve4 = new QwtPlotCurve("Derivative of fit");
    m_curve5 = new QwtPlotCurve("Baseline");
    m_curve6 = new QwtPlotCurve("Control polygon of derivative");
    m_curve7 = new QwtPlotCurve("Zeros of derivative");
    m_curve1->setRenderHint(QwtPlotItem::RenderAntialiased);
    auto sym = new QwtSymbol(QwtSymbol::XCross);
    m_curve1->attach(m_plot);
    m_curve2->setRenderHint(QwtPlotItem::RenderAntialiased);
    sym->setSize(QSize(20,20));
    m_curve2->setSymbol(sym);
    m_curve2->attach(m_plot);
    m_curve3->setRenderHint(QwtPlotItem::RenderAntialiased);
    m_curve3->attach(m_plot);
    m_curve4->attach(m_plot);
    m_curve4->setRenderHint(QwtPlotItem::RenderAntialiased);
    m_curve5->attach(m_plot);
    m_curve5->setRenderHint(QwtPlotItem::RenderAntialiased);
    m_curve6->attach(m_plot);
    m_curve6->setRenderHint(QwtPlotItem::RenderAntialiased);
    m_curve7->attach(m_plot);   
    m_curve7->setRenderHint(QwtPlotItem::RenderAntialiased);
    m_curve7->setStyle(QwtPlotCurve::NoCurve);
    auto sym7 = new QwtSymbol(QwtSymbol::Diamond);
    sym7->setSize(QSize(20,20));
    m_curve7->setSymbol(sym7);

    loadCsvData(csvFile);
    updateSplineApprox();
    mainWidget->show();

    // Create the separate tweaking window
}
void MyMainWindow::loadCsvData(const std::string& csvFile) {
    loadTimeSignal(csvFile, m_dataXs, m_dataYs);
    int numSamples = m_dataXs.size();

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

    m_curve5->setSamples(std::vector<double>({m_dataXs[0], m_dataXs[numSamples-1]}).data(),
                         std::vector<double>({0.0, 0.0}).data(),
                         2);


    // Render the spline approximation.
    std::vector<float> xs;
    std::vector<float> ys;
    bspline_storve::linspace(m_dataXs[0], m_dataXs[numSamples-1], 2000, xs);
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
    
    // Compute spline coefficients of the derivative of the spline fit
    int derDegree = m_degree-1;
    std::vector<float> derCoeffs = bspline_storve::computeDerivativeCoeffs(m_degree, controlPoints, knots);
    for (int i = 0; i < derCoeffs.size(); i++) {
        derCoeffs[i] *= m_derivativeScale;
    }

    // Render the spline derivative
    std::vector<float> ysDer = bspline_storve::renderSpline(derDegree, knots, derCoeffs, xs);
    auto ysDerDouble = stdVectorToDouble(ysDer);
    m_curve4->setSamples(xsDouble.data(), ysDerDouble.data(), xsDouble.size());
    m_curve4->setPen(QPen(QColor(0, 80, 0)));

    // Derivative's control polygon
    std::vector<float> dcpTimes = bspline_storve::controlPolygon(derDegree, knots);
    auto derCoeffsDouble = stdVectorToDouble(derCoeffs);
    auto dcpTimesDouble = stdVectorToDouble(dcpTimes);    
    m_curve6->setSamples(dcpTimesDouble.data(), derCoeffsDouble.data(), cpTimesDouble.size());

    m_curve5->setVisible(m_showBaseline);       
    m_curve1->setVisible(m_showSplineFit);       
    m_curve6->setVisible(m_showDerivPolygon);       
    m_curve4->setVisible(m_showDerivative);

    // Find zeros of the derivative
    int numIterations = 5;
    auto zeroTimes = bspline_storve::computeZeros(derDegree, knots,derCoeffs, numIterations);
    std::cout << "Found " << zeroTimes.size() << " zeros of the derivative\n";
    if (m_showDerivZeros) {
        std::vector<double> derZeroTimes;
        std::vector<double> derZeroValues;
        for (int i = 0; i < zeroTimes.size(); i++) {
            derZeroTimes.push_back(zeroTimes[i]);
            //derZeroValues.push_back(0.00f);
            std::vector<float> tempTimes = {zeroTimes[i]};
            auto tempValues = bspline_storve::renderSpline(m_degree, knots, controlPoints, tempTimes);
            derZeroValues.push_back(tempValues[0]);
        }
        m_curve7->setSamples(derZeroTimes.data(), derZeroValues.data(), zeroTimes.size());
        m_curve7->setVisible(m_showDerivZeros);
    }
    m_plot->replot();
}