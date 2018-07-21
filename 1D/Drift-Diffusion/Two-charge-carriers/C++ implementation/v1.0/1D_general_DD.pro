TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    bernoulli.cpp \
    photogeneration.cpp \
    recombination.cpp \
    Set_diagonals.cpp \
    set_rhs.cpp \
    thomas_tridiag_solve.cpp \
    main.cpp

HEADERS += \
    bernoulli.h \
    parameters.h \
    photogeneration.h \
    recombination.h \
    Set_diagonals.h \
    set_rhs.h \
    thomas_tridiag_solve.h
