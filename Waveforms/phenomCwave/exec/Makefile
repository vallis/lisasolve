CC= g++
C= gcc
CPPFLAGS= -O3 -g -Wall -static 
INCLUDES= -I../include  -I/usr/include -I/home/stas/tools/include
CFLAGS=$(CPPFLAGS)
LIBS= -L/home/stas/tools/lib -L/usr/lib
SDir= ../src/
LDFLAGS= -lm -lfftw3f -lfftw3 -lgsl -lgslcblas 
SOURCES= $(SDir)BBHTemplate.cc  $(SDir)Constants.cc  $(SDir)PhenomCwave.cc\
$(SDir)phenomwf.cc $(SDir)RandomGen.cc \
$(SDir)IndexSet.cc  $(SDir)LISAmotion.cc  $(SDir)LinAlg.cc $(SDir)Matrix.cc $(SDir)fftw++.cc\
$(SDir)ComputeTDI_Freq.cc $(SDir)OrbitalMotion.cc $(SDir)ComputeFisher.cc $(SDir)NoiseModels.cc\
$(SDir)EMRItemplate.cc $(SDir)FastFT.cc $(SDir)FreqAK_RA.cc $(SDir)Kallianpur_Striebel.cc \
$(SDir)MCintegration.cc $(SDir)GSL_vegas.cc

OBJECTS= $(SOURCES:.cc=.o)
EXEC= ComputePhenWave CompTDItry EvalErrorsTest EvalSNRs EvalErrors CompOrb TestFreqAK ComputeEMRIsnr SNR_MBH_SA RS TestMC

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


EvalErrorsTest: $(OBJECTS) EvalErrorsTest.o 
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) EvalErrorsTest.o -o EvalErrorsTest $(LDFLAGS)


#EvalSNRs: $(OBJECTS) EvaluateSNRs.o 
#	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) EvaluateSNRs.o -o EvalSNRs $(LDFLAGS)

EvalSNRs: $(OBJECTS) EvaluateSNRs_miniLISA.o 
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) EvaluateSNRs_miniLISA.o -o EvalSNRs $(LDFLAGS)

SNR_MBH_SA: $(OBJECTS) AverSNRs_miniLISA.o 
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) AverSNRs_miniLISA.o -o SNR_MBH_SA $(LDFLAGS)

EvalErrors: $(OBJECTS) EvalErrorsPhenomC.o 
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) EvalErrorsPhenomC.o -o EvalErrors $(LDFLAGS)

CompOrb: $(OBJECTS) ComputeLISAorbit.o 
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) ComputeLISAorbit.o -o CompOrb $(LDFLAGS)


TestFreqAK: $(OBJECTS) TestAK_inFreq.o 
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) TestAK_inFreq.o -o TestFreqAK $(LDFLAGS)


ComputeEMRIsnr: $(OBJECTS) ComputeEMRIsnrs.o 
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) ComputeEMRIsnrs.o -o ComputeEMRIsnr $(LDFLAGS)


RS: $(OBJECTS) ReconstructSig.o 
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) ReconstructSig.o -o RS $(LDFLAGS)

TestMC: $(OBJECTS) TestMC.o 
	$(CC) $(CPPFLAGS) $(LIBS) $(OBJECTS) TestMC.o -o TestMC $(LDFLAGS)


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

