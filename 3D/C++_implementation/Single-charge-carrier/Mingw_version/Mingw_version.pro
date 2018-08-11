TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -fopenmp

INCLUDEPATH += C:/Eigen
DEPENDPATH += C:/Eigen

SOURCES += \
    ../continuity_p.cpp \
    ../main.cpp \
    ../parameters.cpp \
    ../poisson.cpp \
    ../Utilities.cpp

HEADERS += \
    ../constants.h \
    ../continuity_p.h \
    ../parameters.h \
    ../poisson.h \
    ../Utilities.h

LIBS += -L C:/mingw/ -llibgomp  #need to link this library for it to work with gcc
