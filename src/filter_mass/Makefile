LIBNAME := filter_mass
EXEC := $(LIBNAME) 

FILTLIB := lib$(LIBNAME).a

LIBSRC := filter_mass.c \
		  tree_io.c \
		  reionization.c \
		  save.c

SRC := main.c $(LIBSRC)
OBJS := $(SRC:.c=.o)

LIBOBJS := $(LIBSRC:.c=.o)
FILTLIB := lib$(LIBNAME).a

INCL  := ./Makefile \
		 filter_mass.h \
		 tree_io.h \
		 reionization.h \
		 save.h

ON_CI ?= false
ifeq ($(CI), true)
    ON_CI ?= true
endif

USE-MPI ?= true
ifeq ($(USE-MPI), true)
ifeq ($(ON_CI), false) #  Don't build with MPI if we're on a continuous integration service. 
    OPTS += -DMPI  #  This creates an MPI version that can be used to process files in parallel
    CC = mpicc  # sets the C-compiler
endif
else
    CC = gcc # sets the C-compiler
endif

GSL_DIR := $(shell gsl-config --prefix)
GSL_INCL := $(shell gsl-config --cflags)
GSL_LIBS := $(shell gsl-config --libs)
GSL_LIBDIR := $(GSL_DIR)/lib

OPTIMIZE = -g -O3 -Wall -Werror # optimization and warning flags
#OPTS = -DDEBUG

CFLAGS = $(OPTIMIZE) $(GSL_INCL) $(OPTS) 
LIBS  += -g -lm  $(GSL_LIBS) -lgsl -lgslcblas


all: $(EXEC)

$(EXEC): $(OBJS) $(FILTLIB)
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)

$(FILTLIB): $(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

$(OBJS): $(INCL)

.phony: clean celan celna clena tests
celan celna clena: clean
clean:
	rm -f $(OBJS) $(EXEC) $(FILTLIB)


