TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS +=  -DMKL_LP64 \  #I think ENABLES MKL TO BE USED
                   -openmp  #this works!, b/c when commented out the -openmp flag, a thread call to eigen givess 1 instead of 8 when include this flag.


#QMAKE_CXXFLAGS += -Ox   #is full optimization in Msvc
#makes no speed difference!

INCLUDEPATH += C:/Eigen \
            C:/IntelSWTools/compilers_and_libraries_2018.3.210/windows/mkl/include
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

LIBS += -L"C:/IntelSWTools/compilers_and_libraries_2018.3.210/windows/mkl/lib/intel64_win" -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  #NOTE: there must be no empty  spaces between the -L and the path string!
