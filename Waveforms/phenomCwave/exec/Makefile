CC= g++
C= gcc
CPPFLAGS= -O3 -g -Wall  
INCLUDES= -I../include -I/Users/stba/tools/include 
CFLAGS=$(CPPFLAGS)
LIBS= -L/Users/stba/tools/lib
SDir= ../src/
LDFLAGS= -lm -lfftw3f -lfftw3 -lgsl -lgslcblas 
SOURCES= $(SDir)BBHTemplate.cc  $(SDir)Constants.cc  $(SDir)PhenomCwave.cc\
$(SDir)phenomwf.cc $(SDir)RandomGen.cc \
$(SDir)IndexSet.cc  $(SDir)LISAmotion.cc  $(SDir)LinAlg.cc $(SDir)Matrix.cc $(SDir)fftw++.cc\
$(SDir)ComputeTDI_Freq.cc $(SDir)OrbitalMotion.cc

OBJECTS= $(SOURCES:.cc=.o)
EXEC= ComputePhenWave CompTDItry

.c.o:
	$(C) -c $(INCLUDES) $< -o $(<D)/$(@F)

.cpp.o:
	$(CC) $(CPPFLAGS) $(INCLUDES) -c $(INCLUDES) $< -o $(<D)/$(@F)

.cc.o:
	$(CC) $(CPPFLAGS) -c $(INCLUDES) $< -o $(<D)/$(@F)

all: $(EXEC) ComputePhenWave CompTDItry
	
ComputePhenWave: $(OBJECTS) CompWaveTry.o 
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) CompWaveTry.o -o ComputePhenWave $(LDFLAGS)


CompTDItry: $(OBJECTS) CompTDI_try.o 
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) CompTDI_try.o -o CompTDItry $(LDFLAGS)



#all: SearchX

#FFTw.o: ../src/fftw++.cc ../include/fftw++.h
#	${CC} ${CFLAGS} ${INCLUDES} -c ../src/fftw++.cc

#SearchX.o: SearchX.cc
#	${CC} ${CFLAGS} ${INCLUDES} -c SearchX.cc


#SearchX: SearchX.o FFTw.o 
#	${CC} SearchX.o fftw++.o ${CFLAGS} ${LDFLAGS} -o SearchX

.PHONY: clean

clean:
	\rm -f *.o *~
	\rm -f ../src/*.o ../src/*~
	\rm SearchX

