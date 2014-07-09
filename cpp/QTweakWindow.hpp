#pragma once
#include <QString>
#include <QVector>
#include <QWidget>
#include <QFileInfo>

// Forward declaration.
class VariableEditorWidget;
class QFormLayout;

class QTweakWindow : public QWidget {
Q_OBJECT
public:
    QTweakWindow(QWidget* parent = 0);
    
    void registerVariable(QString name, float* varPtr, float minValue, float maxValue);
    void registerVariable(QString name, int* varPtr, int minValue, int maxValue);
    void registerVariable(QString name, bool* varPtr);
    void registerVariable(QString name, QString* varPtr);
    void registerVariable(QString name, QFileInfo* varPtr);
    
    // Update all values trough the pointers.
    void pushAll();
    
    // Update the widget with values from the pointers.
    void pullAll();

private:
    QVector<VariableEditorWidget*>  m_editors;
    QFormLayout*                    m_formLayout;

private slots:
    void somethingChanged() {
        emit changesAvailable();
    }

signals:
    void changesAvailable();
};