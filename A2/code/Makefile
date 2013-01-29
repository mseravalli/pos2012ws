####################################################################
#                                                                  #
#             Makefile for FIRE Solver Benchmarks                  #
#                                                                  #
####################################################################

PREP = scorep --user

CC = $(PREP) mpicc
CFLAGS = $(METIS_INC) $(PAPI_INC) -Wall -O3 -std=c99 
LIBS =   $(METIS_LIB) $(PAPI_LIB) -lmetis -lm 

LIBPOS=libpos.a
AR = ar
ARFLAGS = rv

SRCS = initialization.c compute_solution.c finalization.c util_read_files.c util_write_files.c test_functions.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: gccg 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gccg: gccg.c $(LIBPOS)
	$(CC) -o $@ $^ $(LIBS) $(CFLAGS)

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o gccg $(LIBPOS)
