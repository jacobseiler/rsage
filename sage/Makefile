ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

EXEC   = sage 

OBJS   = 	./main.o \
			./core_read_parameter_file.o \
			./core_init.o \
			./core_io_tree.o \
			./core_cool_func.o \
			./core_build_model.o \
			./core_save.o \
			./core_mymalloc.o \
			./core_allvars.o \
			./model_infall.o \
			./model_cooling_heating.o \
			./model_starformation_and_feedback.o \
			./model_disk_instability.o \
			./model_reincorporation.o \
			./model_mergers.o \
			./model_misc.o \
			./io/tree_binary.o \
			./self_consistent/selfcon.o \
			./self_consistent/selfcon_grid.o \
			./self_consistent/grid.o \
			./self_consistent/common.o 

INCL   =	./core_allvars.h  \
			./core_proto.h  \
			./core_simulation.h  \
			./io/tree_binary.h \
			./self_consistent/selfcon.h \
			./self_consistent/selfcon_grid.h \
			./self_consistent/grid.h \
			./self_consistent/common.h \
			./Makefile

OPT := -DROOT_DIR='"${ROOT_DIR}"'

USE-MPI = yes# set this if you want to run in 
ifdef USE-MPI
    OPT += -DMPI  #  This creates an MPI version that can be used to process files in parallel
    CC = mpicc  # sets the C-compiler
else
    CC = gcc  # sets the C-compiler
endif

#OPT += -DDEBUG_SELFCON 

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

LIBS   =   -g -lm  $(GSL_LIBS) -lgsl -lgslcblas 

CFLAGS +=   $(OPTIONS) $(OPT) $(OPTIMIZE) 

default: all

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) ./$(EXEC)

all: $(EXEC) 

celan celna clena claen:clean
