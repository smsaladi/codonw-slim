override cflags  = $(CFLAGS) -g

objects = codon_us.o codons.o open_fil.o commline.o coresp.o

CC=cc
CFLAGS= -O -DBSD
LN=ln -f


all: codonw

codonw: $(objects)
	$(CC) $(CFLAGS)  $(objects) -o codonw -lm

clean:
	\rm -f $(objects) codonw

codon_us.o: codon_us.c codonW.h
	$(CC) -c $(CFLAGS) codon_us.c

codons.o: codons.c codonW.h
	$(CC) -c $(CFLAGS) codons.c

coresp.o: coresp.c codonW.h
	$(CC) -c $(CFLAGS) coresp.c

open_fil.o:    open_fil.c codonW.h
	$(CC) -c $(CFLAGS) open_fil.c

commline.o:    commline.c codonW.h
	$(CC) -c $(CFLAGS) commline.c
