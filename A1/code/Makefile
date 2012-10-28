####################################################################
#                                                                  #
#             Makefile for FIRE Solver Benchmarks                  #
#                                                                  #
####################################################################

CC = icc
CFLAGS = -std=c99 -O2 
LIBS = -lm

# Make sure you have loaded the papi module before uncommenting these and 
# include papi.h in the sources
CFLAGS += $(PAPI_INC)
LIBS += $(PAPI_LIB)

SRCS_GCCG = xread.c xwrite.c gccg.c vol2mesh.c
OBJS_GCCG =  $(addsuffix .o, $(basename $(SRCS_GCCG)))

SRCS_BC = binconv.c 
OBJS_BC =  $(addsuffix .o, $(basename $(SRCS_BC)))

all: gccg binconv 

gccg: $(OBJS_GCCG) 
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

binconv: $(OBJS_BC) 
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf *.o gccg binconv *.vtk *.dat
