TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -openmp

INCLUDEPATH += C:/Eigen
DEPENDPATH += C:/Eigen

SOURCES += main.cpp \
    continuity_p.cpp \
    parameters.cpp \
    poisson.cpp \
    Utilities.cpp

HEADERS += \
    constants.h \
    continuity_p.h \
    parameters.h \
    poisson.h \
    Utilities.h
