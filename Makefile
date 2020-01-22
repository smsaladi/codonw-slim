override cflags = $(CFLAGS) -g

objects = codon_all.o codon_idx.o codon_blk.o codons.o commline.o defaults.o

CC=cc
CFLAGS= -O3 -DBSD

all: codonw

codonw: $(objects)
	\mkdir -p bin
	$(CC) $(CFLAGS) $(objects) -o bin/codonw -lm

clean:
	\rm -f $(objects) bin/codonw

codon_all.o: codon_all.c codonW.h
	$(CC) -c $(CFLAGS) codon_all.c

codon_idx.o: codon_idx.c codonW.h
	$(CC) -c $(CFLAGS) codon_idx.c

codon_blk.o: codon_blk.c codonW.h
	$(CC) -c $(CFLAGS) codon_blk.c

codons.o: codons.c codonW.h kseq.h
	$(CC) -c $(CFLAGS) codons.c

commline.o: commline.c codonW.h
	$(CC) -c $(CFLAGS) commline.c

defaults.o: defaults.c codonW.h
	$(CC) -c $(CFLAGS) defaults.c
