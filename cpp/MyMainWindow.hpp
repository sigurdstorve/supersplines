#include <iostream>
#include <string>
#include <vector>
#include <QMainWindow>

// Forward declaration.
class QwtPlotCurve;
class QwtPlot;
class QSlider;
class QPushButton;
class QTweakWindow;

class MyMainWindow : public QMainWindow {
Q_OBJECT
public:
    MyMainWindow(const std::string& csvFile, QWidget* parent = 0); 
private:

    void loadCsvData(const std::string& csvFile);

private slots:
    void setNumControlPoints(int value) {
        std::cout << "Number of control points: " << value << std::endl;
        m_numControlPoints = value;
        updateSplineApprox();
    }
    
    void toggleControlPolygon() {
        m_showControlPolygon = !m_showControlPolygon;
        updateSplineApprox();
    }

    void refresh();

private:
    void updateSplineApprox();
    
    // The raw data
    std::vector<float>      m_dataXs;
    std::vector<float>      m_dataYs;
    
    int                     m_degree;
    int                     m_numControlPoints;
    
    // The curves
    QwtPlotCurve*           m_curve1;
    QwtPlotCurve*           m_curve2;
    QwtPlotCurve*           m_curve3;
    QwtPlotCurve*           m_curve4;
    QwtPlotCurve*           m_curve5;   // baseline
    QwtPlotCurve*           m_curve6;   // control polygon of derivative
    QwtPlotCurve*           m_curve7;   // zeros of derivative

    QwtPlot*                m_plot;
    
    bool                    m_showControlPolygon;
    
    QTweakWindow*           m_tweakWin;
    bool                    m_showDerivative;
    float                   m_derivativeScale;
    bool                    m_showBaseline;
    bool                    m_showDerivPolygon;
    bool                    m_showSplineFit;
    bool                    m_showDerivZeros;
};
