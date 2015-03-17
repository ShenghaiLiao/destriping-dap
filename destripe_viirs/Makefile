INC=-I. -I../destripe/ -I/usr/include
LIBS=-lfftw3f_omp -lfftw3f -lhdf5 -lz -lpthread -lm -lopencv_core
CC=g++
LD=g++
CFLAGS=-g -fopenmp -O2 $(INC)
LDFLAGS=$(LIBS) -fopenmp
TARG=destripe_viirs
OFILES=\
	destripe_viirs.o\
	resample_viirs.o\

HFILES=\
	readwrite_viirs.c\
	resample_viirs.h\
	../destripe/allocate_2d.c\
	../destripe/compute_binary_matrix.c\
	../destripe/compute_meshgrid_dct.c\
	../destripe/compute_nly_filter.c\
	../destripe/compute_poisson_rec.c\
	../destripe/destripe.h\
	../destripe/destripe_main.c\
	../destripe/destripe_main_padframe.c\
	../destripe/fill_restore_viirs_bowtie.c\
	../destripe/get_nif_ndf.c\

all: $(TARG)

$(TARG): $(OFILES)
	$(LD) -o $(TARG) $(OFILES) $(LDFLAGS)

%.o: %.c $(HFILES)
	$(CC) $(CFLAGS) -c $<

install: $(TARG)
	cp $(TARG) /usr/local/bin/

clean:
	rm -f $(OFILES) $(TARG)
