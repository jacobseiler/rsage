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
			./grid_update.o

INCL   =	./core_allvars.h  \
			./core_proto.h  \
			./core_simulation.h  \
			./Makefile

# USE-MPI = yes  # set this if you want to run in parallel

ifdef USE-MPI
    OPT += -DMPI  #  This creates an MPI version that can be used to process files in parallel
    CC = mpicc  # sets the C-compiler
else
    CC = cc  # sets the C-compiler
endif



GSL_INCL = -I/usr/local/include  # make sure your system knows where GSL_DIR is
GSL_LIBS = -L/usr/local/lib
OPTIMIZE = -g -O0 -Wall  # optimization and warning flags

LIBS   =   -g -lm  $(GSL_LIBS) -lgsl -lgslcblas 

CFLAGS =   $(OPTIONS) $(OPT) $(OPTIMIZE) $(GSL_INCL)

default: all

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) ./$(EXEC)

all:  tidy $(EXEC) clean
