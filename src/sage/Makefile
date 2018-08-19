# set this if you want to run in embarrassingly parallel
USE-MPI ?= true
# USE-HDF5 = yes # set this if you want to read in hdf5 trees (requires hdf5 libraries)

ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
LIBS :=
OPT := -DROOT_DIR='"${ROOT_DIR}"'
CCFLAGS := -DGNU_SOURCE -std=gnu99 -fPIC 

LIBNAME = sage
EXEC := $(LIBNAME)

LIBSRC :=  sage.c core_read_parameter_file.c core_init.c core_io_tree.c core_allvars.c \
           core_cool_func.c core_build_model.c core_save.c core_mymalloc.c \
           model_infall.c model_cooling_heating.c model_starformation_and_feedback.c \
           model_disk_instability.c model_reincorporation.c model_mergers.c model_misc.c \
           io/tree_binary.c \
		   self_consistent/selfcon.c self_consistent/selfcon_grid.c self_consistent/grid.c self_consistent/common.c 

SRC := main.c $(LIBSRC)
OBJS := $(SRC:.c=.o)

LIBOBJS := $(LIBSRC:.c=.o)
SAGELIB := lib$(LIBNAME).a

INCL   =	core_allvars.h core_proto.h  core_simulation.h \
			io/tree_binary.h \
			self_consistent/selfcon.h self_consistent/selfcon_grid.h self_consistent/grid.h self_consistent/common.h \
			./Makefile

UNAME := $(shell uname)
ifeq ($(USE-MPI), true)
    OPT += -DMPI  #  This creates an MPI version that can be used to process files in parallel
    CC := mpicc  # sets the C-compiler
else
    # use clang by default on OSX and gcc on linux
    ifeq ($(UNAME), Darwin)
      CC := clang
    else
      CC := gcc
    endif
endif


# No need to do the path + library checks if
# only attempting to clean the build
DO_CHECKS := 1
CLEAN_CMDS := celan celna clean clena
ifneq ($(filter $(CLEAN_CMDS),$(MAKECMDGOALS)),)
  DO_CHECKS := 0
endif

ifeq ($(DO_CHECKS), 1)
  GIT_VERSION ?= $(shell git describe --abbrev=4 --dirty --always --tags)
  OPT += -DVERSION=\"$(GIT_VERSION)\"

  GSL_FOUND ?= $(shell gsl-config --version 2>/dev/null)
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
    GSL_INCL ?= $(shell gsl-config --cflags)
    GSL_LIBDIR ?= $(shell gsl-config --prefix)/lib
    GSL_LIBS   ?= $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR)
  endif
  CFLAGS += $(GSL_INCL)

  OPTIMIZE = -O2 -march=native -mno-fma
  CFLAGS += -g -Wextra -Werror -Wunused-parameter -Wall -Wshadow # optimization and warning flags
  LIBS   =   -g -lm  $(GSL_LIBS) -lgsl -lgslcblas
  CFLAGS += $(OPTIMIZE) $(OPT)
endif

all: $(EXEC) 

$(EXEC): $(OBJS) $(SAGELIB) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)

$(SAGELIB): $(LIBOBJS)
	ar rcs $@ $(LIBOBJS) 

$(OBJS): $(INCL) 

tests: $(EXEC)
	../tests/test_sage.sh

.phony: clean celan celna clena tests
celan celna clena: clean
clean:
	rm -f $(OBJS) $(EXEC) $(SAGELIB)

