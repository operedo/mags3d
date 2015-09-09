
#CC=gcc
CC=icc
#CCFLAGS= -O3 -DOCTAVE -I../include
CCFLAGS= -O3 -xAVX -vec-report -openmp -DGSLIB -I../include
CCLINKS= -lm -lstdc++  

NVCC=/usr/local/cuda-5.5/bin/nvcc 
NVFLAGS= -O3 -arch=sm_20
NVLIBS= -L/usr/local/cuda/lib64 -lcuda -lcudart 

AR= ar
ARFLAGS= cr
RANLIB= ranlib

MAGS3DLIB= libmags3d.a

MAGS3DLIB_OBJECTS= 	genRandomImage.o \
			coordinates.o \
			costFunction.o \
			costFunctionGSLIB.o \
			costFunctionFast.o \
			genMovingAverageImage.o \
			genMovingAverageImageFast.o \
			nscore.o \
			genTargetVariogram.o \
			modifyWeights.o \
			averageWeights.o \
			qsort.o \
			finiteDifferenceGSLIB.o \
			gamv.o \
			gamvCUDA.o \
			gamvCUDAoptimized.o \
			printGSLIB.o

all: lib exe

lib: $(MAGS3DLIB_OBJECTS)  
	$(AR) $(ARFLAGS) $(MAGS3DLIB) $(MAGS3DLIB_OBJECTS)
	$(RANLIB) $(MAGS3DLIB)	
	cp $(MAGS3DLIB) ../lib

exe: mags3d.o mags3dConditioned.o sa-invmags3d.o sa-cuda-invmags3d.o fd-invmags3d.o
	$(CC) $(CCFLAGS) mags3d.o -o mags3d.exe ../lib/$(MAGS3DLIB) $(CCLINKS)
	$(CC) $(CCFLAGS) mags3dConditioned.o -o mags3dConditioned.exe ../lib/$(MAGS3DLIB) $(CCLINKS)
	$(CC) $(CCFLAGS) sa-invmags3d.o -o sa-invmags3d.exe ../lib/$(MAGS3DLIB) $(CCLINKS)
	$(CC) $(CCFLAGS) sa-cuda-invmags3d.o -o sa-cuda-invmags3d.exe ../lib/$(MAGS3DLIB) $(CCLINKS) $(NVLIBS)
	$(CC) $(CCFLAGS) fd-invmags3d.o -o fd-invmags3d.exe ../lib/$(MAGS3DLIB) $(CCLINKS)
	cp *.exe ../bin

clean :
	rm -f *.o *.a *.exe ../bin/*.exe ../lib/*.a

.SUFFIXES: .o .c .cu

.c.o : ; $(CC) -c $(CCFLAGS) $*.c
.cu.o : ; $(NVCC) -c $(NVFLAGS) $*.cu
