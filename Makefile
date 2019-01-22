EXEC   = bin/rsage

SOURCES := src/main.c \
		   src/reion_redshift.c \
		   src/check_ini_files.c

OBJS := $(SOURCES:.c=.o)
INCL   =	src/main.h \
			src/reion_redshift.h \
			src/check_ini_files.h

# Set if you're using self-consistent reionization.
BUILD_RSAGE = false

# Set this to true if you want to run in MPI.
USE-MPI = true

# Determine if we're on continuous integration.
ON_CI := false
ifeq ($(CI), true)
    ON_CI := true
endif

ifeq ($(TRAVIS), true)
    ON_CI := true
endif

ifeq ($(ON_CI), true) #  Don't build with MPI if we're on a continuous integration service. 
    USE-MPI = false
endif

ifeq ($(USE-MPI), true)
    OPT += -DMPI  #  This creates an MPI version that can be used to process files in parallel
    CC = mpicc  # sets the C-compiler
else
    CC = gcc  # sets the C-compiler
endif

# Find the git version of the repo.
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
OPT += -DVERSION=\"$(GIT_VERSION)\"

# Time to add the static libraries.
SAGE_LIB := -Lsrc/sage/ -lsage
ifeq ($(BUILD_RSAGE), true)
	RSAGE_LIB := -Lgrid-model/ -lcifog -Lsrc/filter_mass/ -lfilter_mass
	OPT += -DRSAGE
else
	RSAGE_LIB :=
endif

# Find GSL and add its paths.
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

FFTW3DIR :=/opt/local/include
FFTW_INCL := -I$(FFTW3DIR)
FFTW3_LIBDIR :=/opt/local/lib
FFTW3_LIBS := -L$(FFTW3_LIBDIR) -lfftw3 -Xlinker -rpath -Xlinker $(FFTW3_LIBDIR)

ifeq ($(USE-MPI), true)
ifeq ($(ON_CI), false) #  Don't build with MPI if we're on a continuous integration service.    
    FFTW3_LIBS += -lfftw3_mpi
endif
endif

CFLAGS += $(GSL_INCL) $(FFTW_INCL)

OPTIMIZE = -g -O3 -Wextra -Wall -Wshadow # optimization and warning flags -Wunused-parameter  -Werror 

LIBS   = -lm  $(GSL_LIBS) -lgsl -lgslcblas $(SAGE_LIB) $(RSAGE_LIB) $(FFTW3_LIBS)

CFLAGS += $(OPTIONS) $(OPT) $(OPTIMIZE)

# ==================================
# Shouldn't need to touch below here
# ==================================

# Dont want to redefine some settings in other Makefiles.
# So export them to ensure they're visible.
export USE-MPI
export $(OPTIMIZE)
export $(ON_CI)
export $(GIT_VERSION)
export $(GSL_FOUND)
export $(GSL_INCL)
export $(GSL_LIBDIR)
export $(GSL_LIBS)

# All flags set up, time to compile targets.
default: all 

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL) 

do_tests: compile_sage compile_cifog compile_filter $(EXEC)
	tests/test_sage.sh

clean:
	rm -f $(OBJS) $(EXEC)

tidy:
	rm -f $(OBJS) $(EXEC) 
	$(MAKE) clean -C src/sage/
	$(MAKE) clean -C grid-model/ 
	$(MAKE) clean -C src/filter_mass/

all: compile_sage compile_cifog compile_filter $(EXEC)

sage: compile_sage $(EXEC)

compile_sage:
	$(MAKE) -C src/sage/

compile_cifog:
	$(MAKE) -C grid-model/

compile_filter:
	$(MAKE) -C src/filter_mass/ 

celan celna clena claen: clean
