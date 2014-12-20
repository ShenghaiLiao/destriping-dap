INC=-I. -I../destripe/ -I/usr/include
LIBS=-lfftw3f_omp -lfftw3f -lhdf5 -lsz -lz -lpthread -lm
LDFLAGS=$(LIBS)

# uncomment for ICC
#CC=icc
#LD=icc
#CFLAGS=-openmp -O3 -Zp16 -falign-functions=16 $(INC)
#LDFLAGS += -liomp5

# uncomment for GCC
CC=gcc
LD=gcc
CFLAGS=-g -Wall -fopenmp -O2 $(INC)
LDFLAGS += -fopenmp


PROG=destripe_viirs

OFILES=\
	destripe_viirs.o\

HFILES=\
	readwrite_viirs.c\
	resample_viirs.h\
	../destripe/destripe.h\
	../destripe/allocate_2d.c\
	../destripe/compute_binary_matrix.c\
	../destripe/compute_meshgrid_dct.c\
	../destripe/compute_nly_filter.c\
	../destripe/compute_poisson_rec.c\
	../destripe/destripe_main.c\
	../destripe/get_composite_number.c\

all: $(PROG)

$(PROG): $(OFILES)
	$(LD) -o $(PROG) $(OFILES) $(LDFLAGS)


%.o: %.c $(HFILES)
	$(CC) $(CFLAGS) -c $<

install: $(PROG)
	cp $(PROG) /usr/local/bin/

clean:
	rm -f $(OFILES) $(PROG)
