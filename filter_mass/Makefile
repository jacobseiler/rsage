EXECS  := create_filter_mass 

OBJS  := main.o \
		 read_parameter_file.o \
		 tree_io.o \
		 reionization.o \
		 save.o
 	 
EXECS_OBJS := create_filter_mass.o 

INCL  := ./Makefile \
		 main.h \
		 read_parmater_file.h \
		 tree_io.h \
		 reionization.h \
		 save.h

USE-MPI = yes

ifdef USE-MPI
    OPTS = -DMPI  #  This creates an MPI version that can be used to process files in parallel
    CC = mpicc  # sets the C-compiler
else
    CC = gcc  # sets the C-compiler
endif

ECHO = /bin/echo

GSL_DIR := $(shell gsl-config --prefix)
GSL_INCL := $(shell gsl-config --cflags)
GSL_LIBS := $(shell gsl-config --libs)
GSL_LIBDIR := $(GSL_DIR)/lib

OPTIMIZE = -g -O0 -Wall -Werror # optimization and warning flags
#OPTS = -DDEBUG

CFLAGS = $(OPTIMIZE) $(GSL_INCL) $(OPTS) 
LIBS  += -g -lm  $(GSL_LIBS) -lgsl -lgslcblas

all: $(EXECS)

create_filter_mass: $(OBJS) main.o
	$(CC) $(CFLAGS) $^ $(LIBS) -Xlinker -rpath -Xlinker $(GSL_LIBDIR) -o  $@

clean:
	rm -f $(OBJS) $(EXECS_OBJS) $(EXECS)

.PHONY: all clean clena celan celna

celan celna clena claen:clean
