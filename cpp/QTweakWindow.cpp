#include <stdexcept>
#include <QFormLayout>
#include "QTweakWindow.hpp"
#include "VariableEditorWidget.hpp"

QTweakWindow::QTweakWindow(QWidget* parent) :
        QWidget(parent) {
    //layout = new QVBoxLayout;
    m_formLayout = new QFormLayout;
    setLayout(m_formLayout);
}        

void QTweakWindow::registerVariable(QString name, float* varPtr, float minValue, float maxValue) {
    VariableEditorWidget* varEditor = new FloatEditorWidget(varPtr, minValue, maxValue);
    m_formLayout->addRow(name, varEditor);
    connect(varEditor, SIGNAL(valueHasChanged()), this, SLOT(somethingChanged()));
    m_editors.append(varEditor);
}
void QTweakWindow::registerVariable(QString name, int* varPtr, int minValue, int maxValue) {
    VariableEditorWidget* varEditor = new IntEditorWidget(varPtr, minValue, maxValue);
    m_formLayout->addRow(name, varEditor);
    connect(varEditor, SIGNAL(valueHasChanged()), this, SLOT(somethingChanged()));
    m_editors.append(varEditor);
}
void QTweakWindow::registerVariable(QString name, bool* varPtr) {
    VariableEditorWidget* varEditor = new BoolEditorWidget(varPtr);
    m_formLayout->addRow(name, varEditor);
    connect(varEditor, SIGNAL(valueHasChanged()), this, SLOT(somethingChanged()));
    m_editors.append(varEditor);
}
void QTweakWindow::registerVariable(QString name, QString* varPtr) {
    VariableEditorWidget* varEditor = new QStringEditorWidget(varPtr);
    m_formLayout->addRow(name, varEditor);
    connect(varEditor, SIGNAL(valueHasChanged()), this, SLOT(somethingChanged()));
    m_editors.append(varEditor);
}
void QTweakWindow::registerVariable(QString name, QFileInfo* varPtr) {
    VariableEditorWidget* varEditor = new QFileInfoEditorWidget(varPtr);
    m_formLayout->addRow(name, varEditor);
    connect(varEditor, SIGNAL(valueHasChanged()), this, SLOT(somethingChanged()));
    m_editors.append(varEditor);
}

void QTweakWindow::pushAll() {
    for (auto it = m_editors.begin(); it != m_editors.end(); ++it) {
        (*it)->push();
    }
}

void QTweakWindow::pullAll() {
    for (auto it = m_editors.begin(); it != m_editors.end(); ++it) {
        (*it)->pull();
    }
}
