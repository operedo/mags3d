
MAGS3DLIB= libmags3d.a

CC=icc
CCFLAGS= -O3 -openmp -I../include
CCLINKS= ../lib/$(MAGS3DLIB) -lm
NVLIBS= -L/usr/local/cuda/lib64 -lcuda -lcudart 

all: 	testGamv
#	testGenRandomImage \
#	testCoordinates \

# not finished
testGenRandomImage: testGenRandomImage.o
	$(CC) $(CCFLAGS) testGenRandomImage.o -o testGenRandomImage.exe $(CCLINKS)

# not finished
testCoordinates: testCoordinates.o
	$(CC) $(CCFLAGS) testCoordinates.o -o testCoordinates.exe $(CCLINKS)

testGamv: testGamv.o
	$(CC) $(CCFLAGS) testGamv.o -o testGamv.exe $(CCLINKS) $(NVLIBS)


.PHONY: clean

clean :
	rm -f *.exe *.o

.SUFFIXES: .o .c

.c.o : ; $(CC) -c $(CCFLAGS) $*.c

