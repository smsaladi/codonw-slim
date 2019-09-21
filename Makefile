override cflags = $(CFLAGS) -g

objects = codon_us.o codons.o commline.o coresp.o defaults.o

CC=cc
CFLAGS= -O3 -DBSD

all: codonw

codonw: $(objects)
	\mkdir -p bin
	$(CC) $(CFLAGS) $(objects) -o bin/codonw -lm

clean:
	\rm -f $(objects) bin/codonw

codon_us.o: codon_us.c codonW.h
	$(CC) -c $(CFLAGS) codon_us.c

codons.o: codons.c codonW.h kseq.h
	$(CC) -c $(CFLAGS) codons.c

coresp.o: coresp.c codonW.h
	$(CC) -c $(CFLAGS) coresp.c

commline.o: commline.c codonW.h
	$(CC) -c $(CFLAGS) commline.c

defaults.o: defaults.c codonW.h
	$(CC) -c $(CFLAGS) defaults.c
