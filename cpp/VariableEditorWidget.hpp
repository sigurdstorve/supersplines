#pragma once
#include <QBoxLayout>
#include <QLabel>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QLineEdit>
#include <QVector>

class VariableEditorWidget : public QWidget {
Q_OBJECT
public:
    VariableEditorWidget(QWidget* parent = 0) :
        QWidget(parent) {
    }
    
    // Push changes through pointer.
    virtual void push()     = 0;
    
    // Pull changes from pointer.
    virtual void pull()     = 0;

signals:
    void valueHasChanged();
};

class FloatEditorWidget : public VariableEditorWidget {
Q_OBJECT
public:
    FloatEditorWidget(float* varPtr, float minValue, float maxValue, QWidget* parent = 0) :
            VariableEditorWidget(parent) {
        QHBoxLayout* layout = new QHBoxLayout;
        setLayout(layout);
        
        m_spinBox = new QDoubleSpinBox;
        m_spinBox->setRange(minValue, maxValue);
        m_spinBox->setSingleStep(1);
        connect(m_spinBox, SIGNAL(valueChanged(double)), this, SLOT(valueChanged(double)));
        m_varPtr = varPtr;
        layout->addWidget(m_spinBox);
        pull();
    }
        
    virtual void push() {
        *m_varPtr = static_cast<float>(m_cachedValue);
    }

    virtual void pull() {
        m_cachedValue = static_cast<double>(*m_varPtr);
        m_spinBox->setValue(m_cachedValue);
    }

private slots:
    void valueChanged(double newValue) {
        m_cachedValue = newValue;
        emit valueHasChanged();
    }

private:
    float*              m_varPtr;
    double              m_cachedValue;
    QDoubleSpinBox*     m_spinBox;
};

class IntEditorWidget : public VariableEditorWidget {
Q_OBJECT
public:
    IntEditorWidget(int* varPtr, int minValue, int maxValue, QWidget* parent = 0) :
            VariableEditorWidget(parent) {
        QHBoxLayout* layout = new QHBoxLayout;
        setLayout(layout);
        
        m_spinBox = new QSpinBox;
        m_spinBox->setRange(minValue, maxValue);
        m_spinBox->setSingleStep(1);
        connect(m_spinBox, SIGNAL(valueChanged(int)), this, SLOT(valueChanged(int)));
        m_varPtr = varPtr;
        layout->addWidget(m_spinBox);
        pull();
    }
        
    virtual void push() {
        *m_varPtr = m_cachedValue;
    }

    virtual void pull() {
        m_cachedValue = *m_varPtr;
        m_spinBox->setValue(m_cachedValue);
    }

private slots:
    void valueChanged(int newValue) {
        m_cachedValue = newValue;
        emit valueHasChanged();
    }

private:
    int*                m_varPtr;
    int                 m_cachedValue;
    QSpinBox*           m_spinBox;
};

class BoolEditorWidget : public VariableEditorWidget {
Q_OBJECT
public:
    BoolEditorWidget(bool* varPtr, QWidget* parent = 0) :
            VariableEditorWidget(parent) {
        QHBoxLayout* layout = new QHBoxLayout;
        setLayout(layout);
        
        m_checkBox = new QCheckBox;
        connect(m_checkBox, SIGNAL(stateChanged(int)), this, SLOT(valueChanged(int)));
        m_varPtr = varPtr;
        layout->addWidget(m_checkBox);
        pull();
    }
        
    virtual void push() {
        *m_varPtr = m_cachedValue;
    }

    virtual void pull() {
        m_cachedValue = *m_varPtr;
        m_checkBox->setChecked(m_cachedValue);
    }

private slots:
    void valueChanged(int newValue) {
        m_cachedValue = (newValue == 2);
        emit valueHasChanged();
    }

private:
    bool*               m_varPtr;
    bool                m_cachedValue;
    QCheckBox*          m_checkBox;
};

class QStringEditorWidget : public VariableEditorWidget {
Q_OBJECT
public:
    QStringEditorWidget(QString* varPtr, QWidget* parent = 0) :
            VariableEditorWidget(parent) {
        QHBoxLayout* layout = new QHBoxLayout;
        setLayout(layout);
        
        m_lineEdit = new QLineEdit;
        connect(m_lineEdit, SIGNAL(textChanged(QString)), this, SLOT(valueChanged(QString)));
        m_varPtr = varPtr;
        layout->addWidget(m_lineEdit);
        pull();
    }
        
    virtual void push() {
        *m_varPtr = m_cachedValue;
    }

    virtual void pull() {
        m_cachedValue = *m_varPtr;
        m_lineEdit->setText(m_cachedValue);
    }

private slots:
    void valueChanged(QString newValue) {
        m_cachedValue = newValue;
        emit valueHasChanged();
    }

private:
    QString*            m_varPtr;
    QString             m_cachedValue;
    QLineEdit*          m_lineEdit;
};
