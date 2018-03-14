EXECS  := create_filter_mass 

OBJS  := filter_mass/main.o \
		 filter_mass/read_parameter_file.o \
		 filter_mass/tree_io.o \
		 filter_mass/reionization.o \
		 filter_mass/save.o
 	 
EXECS_OBJS := filter_mass/create_filter_mass.o 

INCL  := ./Makefile \
		 filter_mass/main.h \
		 filter_mass/read_parmater_file.h \
		 filter_mass/tree_io.h \
		 filter_mass/reionization.h \
		 filter_mass/save.h

ECHO = /bin/echo

CC := cc  # sets the C-compiler

GSL_DIR := $(shell gsl-config --prefix)
GSL_INCL := $(shell gsl-config --cflags)
GSL_LIBS := $(shell gsl-config --libs)
GSL_LIBDIR := $(GSL_DIR)/lib

OPTIMIZE = -g -O0 -Wall -Werror # optimization and warning flags
#OPTS = -DDEBUG

CFLAGS = $(OPTIMIZE) $(GSL_INCL) $(OPTS) 
LIBS  += -g -lm  $(GSL_LIBS) -lgsl -lgslcblas

all: $(EXECS)

create_filter_mass: $(OBJS) filter_mass/main.o
	$(CC) $(CFLAGS) $^ $(LIBS) -Xlinker -rpath -Xlinker $(GSL_LIBDIR) -o  $@

clean:
	rm -f $(OBJS) $(EXECS_OBJS) $(EXECS)

.PHONY: all clean clena celan celna

celan celna clena claen:clean
