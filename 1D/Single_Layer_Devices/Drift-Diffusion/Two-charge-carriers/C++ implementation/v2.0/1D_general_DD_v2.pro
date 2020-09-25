TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
#QMAKE_CXXFLAGS_RELEASE += -Ox  //Ox is "full optimization" for Msvc, seems no difference in speedfrom the default -O2

SOURCES += \
    main.cpp \
    photogeneration.cpp \
    recombination.cpp \
    thomas_tridiag_solve.cpp \
    poisson.cpp \
    continuity_n.cpp \
    continuity_p.cpp \
    parameters.cpp \
    Utilities.cpp


HEADERS += \
    photogeneration.h \
    recombination.h \
    thomas_tridiag_solve.h \
    poisson.h \
    continuity_n.h \
    continuity_p.h \
    constants.h \
    parameters.h \
    Utilities.h
