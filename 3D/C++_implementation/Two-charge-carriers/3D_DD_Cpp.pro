TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -openmp

INCLUDEPATH += C:/Eigen
DEPENDPATH += C:/Eigen

SOURCES += main.cpp \
    continuity_n.cpp \
    continuity_p.cpp \
    parameters.cpp \
    photogeneration.cpp \
    poisson.cpp \
    recombination.cpp \
    Utilities.cpp

HEADERS += \
    constants.h \
    continuity_n.h \
    continuity_p.h \
    parameters.h \
    photogeneration.h \
    poisson.h \
    recombination.h \
    Utilities.h
