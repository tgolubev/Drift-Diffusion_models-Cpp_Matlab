TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    bernoulli.cpp \
    photogeneration.cpp \
    recombination.cpp \
    thomas_tridiag_solve.cpp \
    main.cpp \
    poisson.cpp \
    continuity_n.cpp \
    continuity_p.cpp \
    parameters.cpp

HEADERS += \
    bernoulli.h \
    photogeneration.h \
    recombination.h \
    thomas_tridiag_solve.h \
    poisson.h \
    continuity_n.h \
    continuity_p.h \
    constants.h \
    parameters.h
