# Comment lines
# Here we define compiler option, libraries and the target
FLAGS = -O2 -Wall

# Here we make the executable file
SRCS =   main.cpp bernoulli.cpp photogeneration.cpp recombination.cpp Set_diagonals.cpp set_rhs.cpp thomas_tridiag_solve.cpp
OBJS = $(subst .cpp,.o,$(SRCS))
all: DD

# Whereas here we create the object file
DD: $(OBJS) parameters.h #need to put this explicitely b/ doesn't have an associated .cpp
	g++ ${FLAGS} -o DD $(OBJS) 
#	export OMP_NUM_THREADS=20
#	./md

%.o: %.cpp  parameters.h
	g++ $(FLAGS) -c $<


# Clean
clean:
	rm *.o ./DD
