LIB := libfilter_mass.a 
EXEC := filter_mass 

OBJS  := filter_mass.o \
		 tree_io.o \
		 reionization.o \
		 save.o
 	 
EXECS_OBJS := create_filter_mass.o 

INCL  := ./Makefile \
		 filter_mass.h \
		 tree_io.h \
		 reionization.h \
		 save.h

ON_CI ?= false
ifeq ($(CI), true)
    ON_CI ?= true
endif

USE-MPI ?= yes
ifdef USE-MPI
ifeq ($(ON_CI), false) #  Don't build with MPI if we're on a continuous integration service. 
    OPT += -DMPI  #  This creates an MPI version that can be used to process files in parallel
    CC = mpicc  # sets the C-compiler
endif
else
    CC = gcc # sets the C-compiler
endif

GSL_DIR := $(shell gsl-config --prefix)
GSL_INCL := $(shell gsl-config --cflags)
GSL_LIBS := $(shell gsl-config --libs)
GSL_LIBDIR := $(GSL_DIR)/lib

OPTIMIZE = -g -O0 -Wall -Werror # optimization and warning flags
#OPTS = -DDEBUG

CFLAGS = $(OPTIMIZE) $(GSL_INCL) $(OPTS) 
LIBS  += -g -lm  $(GSL_LIBS) -lgsl -lgslcblas

default: lib

$(LIB) : $(OBJS) 
	$(CC) -c $(CFLAGS) $(LIBS) -Xlinker -rpath -Xlinker $(GSL_LIBDIR) 
	ar rcs $(LIB) $(OBJS)

lib: $(LIB)

clean:
	rm -f $(OBJS) $(EXECS_OBJS) $(EXECS) $(LIB)

.PHONY: all clean clena celan celna

celan celna clena claen:clean
