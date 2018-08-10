ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

EXEC   = rsage 

OBJS   = 	./main.o 

INCL   =	./main.h \
			./Makefile

BUILD_SAGE = yes 
BUILD_RSAGE = yes 

OPT := -DROOT_DIR='"${ROOT_DIR}"'

ON_CI := false
ifeq ($(CI), true)
    ON_CI := true
endif

ifeq ($(TRAVIS), true)
    ON_CI := true
endif

USE-MPI = yes# set this if you want to run in 

ifdef USE-MPI
ifeq ($(ON_CI), false) #  Don't build with MPI if we're on a continuous integration service. 
    OPT += -DMPI  #  This creates an MPI version that can be used to process files in parallel
    CC = mpicc  # sets the C-compiler
endif
else
    CC = gcc  # sets the C-compiler
endif

GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
OPT += -DVERSION=\"$(GIT_VERSION)\"

# Adds the SAGE library. 
ifdef BUILD_SAGE
	SAGE_LIB := -Lsage/ -lsage
else
	SAGE_LIB :=
endif

# Adds the cifog library. 
ifdef BUILD_RSAGE
	RSAGE_LIB := -L../grid-model/ -lcifog -Lfilter_mass/ -lfilter_mass
else
	RSAGE_LIB :=
endif

GSL_FOUND := $(shell gsl-config --version 2>/dev/null)
ifndef GSL_FOUND
  $(warning GSL not found in path - please install GSL before installing SAGE (or, update the PATH environment variable such that "gsl-config" is found))
  $(warning Assuming GSL *might* be in $(GSL_DIR) and trying to compile)
  # if the automatic detection fails, set GSL_DIR appropriately
  GSL_DIR := /opt/local
  GSL_INCL := -I$(GSL_DIR)/include  
  GSL_LIBDIR := $(GSL_DIR)/lib
  # since GSL is not in PATH, the runtime environment might not be setup correctly either
  # therefore, adding the compiletime library path is even more important (the -Xlinker bit)
  GSL_LIBS := -L$(GSL_LIBDIR) -lgsl -lgslcblas -Xlinker -rpath -Xlinker $(GSL_LIBDIR) 
else
  # GSL is probably configured correctly, pick up the locations automatically
  GSL_INCL := $(shell gsl-config --cflags)
  GSL_LIBDIR := $(shell gsl-config --prefix)/lib
  GSL_LIBS   := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR)
endif
CFLAGS += $(GSL_INCL)

OPTIMIZE = -g -O3 -Wextra -Werror -Wunused-parameter -Wall -Wshadow # optimization and warning flags

LIBS   =   -g -lm  $(GSL_LIBS) -lgsl -lgslcblas $(SAGE_LIB) $(RSAGE_LIB) 

CFLAGS += $(OPTIONS) $(OPT) $(OPTIMIZE) 

export USE-MPI

default: all 

$(EXEC): $(OBJS)
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)

$(OBJS): $(INCL) 

tests: $(EXEC)
	../tests/test_sage.sh

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) ./$(EXEC)

all: compile_sage compile_cifog compile_filter $(EXEC)

sage: compile_sage $(EXEC)

compile_sage:
	$(MAKE) -C sage/

compile_cifog:
	$(MAKE) -C ../grid-model/

compile_filter:
	$(MAKE) -C filter_mass/ 

celan celna clena claen: clean
