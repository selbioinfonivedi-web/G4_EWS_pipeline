# compiler 
#CC        = /opt/ibmcmp/vac/6.0/bin/xlc
CC        = gcc
CFLAGS  += -O3  -Wall

# target macros
PhiOBJ	= normal.o stats.o maxChi.o main.o queue.o graphCode.o fasta.o phylip.o mem.o misc.o pairScore.o seqManip.o global.o
ProfOBJ = normal.o stats.o fasta.o phylip.o mem.o misc.o pairScore.o seqManip.o global.o profile.o

# targets
default : Phi Profile

# implicit construction rule
%.o : %.c
	$(CC) -c  $(CFLAGS) $< -o $@


Phi: $(PhiOBJ) 
	$(CC) $(CFLAGS)  -o Phi $(PhiOBJ)   -lm $(LDFLAGS)
	cp Phi ../
	cd ppma_2_bmp && ${MAKE}
	cp ppma_2_bmp/ppma_2_bmp ../


Profile: $(ProfOBJ)
	$(CC) $(CFLAGS) -o Profile $(ProfOBJ) -lm $(MYLIB) $(LDFLAGS)
	cp Profile ../
clean : FORCE
	rm -f *.o
	rm -f Phi
	rm -f Profile
FORCE:



