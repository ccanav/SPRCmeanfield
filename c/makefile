radau5.o: bvia.h radau5.f
	gfortran -L/usr/local/g2c/lib   -c radau5.f
decsol.o: bvia.h decsol.f
	gfortran  -L/usr/local/g2c/lib -c decsol.f
pviaprcout.o : bvia.h pviaprcout.c
	gcc  -c pviaprcout.c 
aprcvia.o : avia.h aprcvia.c
	gcc  -c aprcvia.c 
aviaprcout.o : avia.h aviaprcout.c
	gcc  -c aviaprcout.c 
aviaprcdriver.o : avia.h aprcvia.c aviaprcout.c aviaprcdriver.c
	gcc  -c aviaprcdriver.c 
aprcvia:  radau5.o decsol.o aviaprcout.c aviaprcdriver.c aprcvia.c avia.h cblock.h aviaprcout.o aviaprcdriver.o 
	gfortran   -o aprcvia radau5.o decsol.o  aviaprcout.o  aviaprcdriver.o -lm 
mprcvia.o : bvia.h mprcvia.c
	gcc  -c mprcvia.c 
mviaprcdriver.o : bvia.h mprcvia.c pviaprcout.c mviaprcdriver.c
	gcc  -c mviaprcdriver.c 
mprcvia:  radau5.o decsol.o pviaprcout.c mviaprcdriver.c mprcvia.c bvia.h cblock.h pviaprcout.o mviaprcdriver.o 
	gfortran   -o mprcvia radau5.o decsol.o  pviaprcout.o  mviaprcdriver.o -lm 
dviaprcdriver.o : bvia.h dprcvia.c pviaprcout.c dviaprcdriver.c
	gcc  -c dviaprcdriver.c 
freqpred.o : bvia.h dprcvia.c pviaprcout.c freqpred.c
	gcc  -c freqpred.c 
dprcvia:  radau5.o decsol.o pviaprcout.c dviaprcdriver.c dprcvia.c bvia.h cblock.h pviaprcout.o dviaprcdriver.o 
	gfortran   -o dprcvia radau5.o decsol.o  pviaprcout.o  dviaprcdriver.o -lm 
sprcvia:  radau5.o decsol.o pviaprcout.c sviaprcdriver.c sprcvia.c bvia.h cblock.h pviaprcout.o sviaprcdriver.o 
	gfortran   -o sprcvia radau5.o decsol.o  pviaprcout.o  sviaprcdriver.o -lm 
freqpred:  radau5.o decsol.o pviaprcout.c freqpred.c dprcvia.c bvia.h cblock.h pviaprcout.o freqpred.o 
	gfortran   -o freqpred radau5.o decsol.o  pviaprcout.o  freqpred.o -lm 
