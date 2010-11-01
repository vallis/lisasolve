%module FastBinary
%{
#include "FastBinary.h"
%}

// the following three lines allow the automatic conversion of a numpy array
// into a pointer to an array of doubles and the long length of the array
// each pair of these that appears in the library header file (e.g., BBHChallenge1.hh)
// needs to be named explicitly

%include numpy_typemaps.i
%typemap(in) (double *XLS,long XLSlen) = (double *numpyarray,long numpyarraysize);
%typemap(in) (double *XSL,long XSLlen) = (double *numpyarray,long numpyarraysize);
%typemap(in) (double *YLS,long YLSlen) = (double *numpyarray,long numpyarraysize);
%typemap(in) (double *YSL,long YSLlen) = (double *numpyarray,long numpyarraysize);
%typemap(in) (double *ZLS,long ZLSlen) = (double *numpyarray,long numpyarraysize);
%typemap(in) (double *ZSL,long ZSLlen) = (double *numpyarray,long numpyarraysize);

// the following include defines all the calls that need to be wrapped by SWIG
// the include at the top serves only to provide definitions for this one;
// since it is between braces %{ %}, its content goes straight into
// the CPP wrapper file

%include "FastBinary.h"

// and here we can add some Python to the interface code

%pythoncode %{
import lisaxml
# import lisaxml2 as lisaxml

import numpy
import FrequencyArray
import math
import time
import sys

lisaxml.SourceClassModules['FastGalacticBinary'] = 'FastBinary'

# define the interface class

class FastGalacticBinary(lisaxml.Source):
    outputlist = ( ('Frequency',                  'Hertz',     None, 'GW frequency'),
                   ('FrequencyDerivative',        'Hertz/s',   0,    'GW frequency derivative'),
                   ('EclipticLatitude',           'Radian',    None, 'standard ecliptic latitude'),
                   ('EclipticLongitude',          'Radian',    None, 'standard ecliptic longitude'),
                   ('Amplitude',                  '1',         None, 'dimensionless GW amplitude'),
                   ('Inclination',                'Radian',    None, 'standard source inclination'),
                   ('Polarization',               'Radian',    None, 'standard source polarization'),
                   ('InitialPhase',               'Radian',    None, 'GW phase at t = 0') )
    
    FastBinaryCache = {}
    
    # TO DO: the initialization may be a useful extension of lisaxml.Source
    def __init__(self,name='',init=None):
        super(FastGalacticBinary, self).__init__('FastGalacticBinary',name)
        
        # try to initialize the parameters from "init"
        if isinstance(init,(tuple,list)) and len(init) == len(self.outputlist):
            for i,value in enumerate(init):
                setattr(self,self.outputlist[i][0],value)
        elif isinstance(init,dict):
            for par in self.outputlist:
                if par[0] in init:
                    setattr(self,par[0],init[par[0]])
                elif par[2] != None:
                    setattr(self,par[0],par[2])
    
    def fourier(self,simulator='synthlisa',table=None,T=6.2914560e7,dt=15,mpi=False):
        if table == None:
            return self.onefourier(simulator,T=T,dt=dt)
        elif mpi == False:
            NFFT = int(T/dt)
            buf = tuple(FrequencyArray.FrequencyArray(numpy.zeros(NFFT/2+1,dtype=numpy.complex128),kmin=0,df=1.0/T) for i in range(3))
            
            lng = len(table); cnt = 0; time0 = time.time()
            for line in table:
                (self.Frequency,
                 self.FrequencyDerivative,
                 self.EclipticLatitude,
                 self.EclipticLongitude,
                 self.Amplitude,
                 self.Inclination,
                 self.Polarization,
                 self.InitialPhase) = line[:]
                
                self.onefourier(simulator,buffer=buf,T=T,dt=dt)
                
                cnt = cnt + 1
                if cnt % 10000 == 0:
                    time1 = time.time()
                    print "\r%d/%d sources (%d/s), ETA %d s                    " % (cnt,lng,cnt/(time1-time0),(lng-cnt)/(cnt/(time1-time0))),
                    sys.stdout.flush()
            time1 = time.time()
            print "\r%d finished, %d s elapsed (%d/s)" % (cnt,time1-time0,cnt/(time1-time0))
            
            return buf
        else:
            from mpi4py import MPI
            
            size, rank = MPI.COMM_WORLD.Get_size(), MPI.COMM_WORLD.Get_rank()
            slaves = size - 1
            
            send_bufs = [numpy.zeros(8,'d') for s in range(slaves)]
            send_reqs = [MPI.COMM_WORLD.Send_init((send_bufs[s], MPI.DOUBLE), s+1, 0) for s in range(slaves)]
            
            recv_bufs = [numpy.zeros(1,'d') for s in range(slaves)]
            recv_reqs = [MPI.COMM_WORLD.Recv_init((recv_bufs[s], MPI.DOUBLE), s+1, 1) for s in range(slaves)]
            
            busy = [-1 for s in range(slaves)]
            
            ln = len(table); cnt = 0
            import countdown; c = countdown.countdown(ln,10000)
            while cnt < ln:
                for s in range(slaves):                 
                    if busy[s] == -1 and cnt < ln:
                        send_bufs[s][:] = table[cnt][:]
                        send_reqs[s].Start()
                        recv_reqs[s].Start()
                        
                        # print "Asking slave %d to do source %d" % (s+1,cnt); sys.stdout.flush()
                        
                        busy[s] = cnt
                        cnt = cnt + 1
                
                if sum(x == -1 for x in busy) != len(busy):     # at least one of busy[s] is not -1
                    f = MPI.Prequest.Waitany(recv_reqs)
                    
                    # print "Completed source %d from slave %d" % (busy[f],f+1); sys.stdout.flush()
                    
                    busy[f] = -1
                    c.status()
            c.end()
            
            for s in range(0,slaves):
                send_bufs[s][0] = 0.0   # terminate by sending a zero as first parameter
                send_reqs[s].Start()    
            
            NFFT = int(T/dt)
            buf = tuple(FrequencyArray.FrequencyArray(numpy.zeros(NFFT/2+1,dtype=numpy.complex128),kmin=0,df=1.0/T) for i in range(3))
            rbuf = numpy.zeros(NFFT/2+1,dtype=numpy.complex128)
            
            for s in range(slaves):
                for i in range(3):
                    MPI.COMM_WORLD.Recv((rbuf,MPI.DOUBLE_COMPLEX),s+1,i)
                    buf[i][:] += rbuf[:]
            
            # something strange going on with mpi4py's Reduce. Will use simple communication.
            # for i in range(3):
            #     MPI.COMM_WORLD.Reduce(None,(buf[i],MPI.DOUBLE_COMPLEX),MPI.SUM,0)
            
            return buf
    
    year = 3.15581498e7
    fstar = 0.00954269032
    
    def slave(self,simulator='synthlisa',T=6.2914560e7,dt=15):
        from mpi4py import MPI
        
        NFFT = int(T/dt)
        buf = tuple(FrequencyArray.FrequencyArray(numpy.zeros(NFFT/2+1,dtype=numpy.complex128),kmin=0,df=1.0/T) for i in range(3))
        
        recv_buf, send_buf = numpy.zeros(8,'d'), numpy.zeros(1,'d')
        recv_req, send_req = MPI.COMM_WORLD.Recv_init((recv_buf, MPI.DOUBLE), 0, 0), MPI.COMM_WORLD.Send_init((send_buf, MPI.DOUBLE), 0, 1)
        
        while True:
            recv_req.Start()
            MPI.Prequest.Wait(recv_req)
            
            if recv_buf[0] == 0.0:
                break
            else:
                self.onefourier(simulator,vector=recv_buf,buffer=buf,T=T,dt=dt)
                send_req.Start()
        
        for i in range(3):
            MPI.COMM_WORLD.Send((buf[i],MPI.DOUBLE_COMPLEX),0,i)
        
        # for i in range(3):
        #     MPI.COMM_WORLD.Reduce((buf[i],MPI.DOUBLE_COMPLEX),None,MPI.SUM,0)
    
    def onefourier(self,simulator='synthlisa',vector=None,buffer=None,T=6.2914560e7,dt=15):
        if vector != None:
            self.Frequency,self.Amplitude = vector[0],vector[4]
        
        if T/self.year <= 1.0:
            mult = 1
        elif T/self.year <= 2.0:
            mult = 2
        elif T/self.year <= 4.0:
            mult = 4
        else:
            mult = 8
        
        if self.Frequency > 0.1:
            N = 1024*mult
        elif self.Frequency > 0.03:
            N = 512*mult
        elif self.Frequency > 0.01:
            N = 256*mult
        elif self.Frequency > 0.001:
            N = 64*mult
        else:
            N = 32*mult            
        
        Sm = AEnoise(self.Frequency) / (4.0 * math.sin(self.Frequency/self.fstar) * math.sin(self.Frequency/self.fstar))
        Acut = self.Amplitude * math.sqrt(T/Sm)
        
        M = int(math.pow(2.0,1 + int(math.log(Acut)/math.log(2.0))));
        
        # corrected per Neil's 2009-07-28 changes to Fast_Response3.c
        M = N = min(8192,max(M,N))
        
        # cache FastResponse objects
        if (N,T,dt) not in FastGalacticBinary.FastBinaryCache:
            fastbin = FastResponse(N,T,dt)
            fastbin.XLS = numpy.zeros(2*M,'d'); fastbin.XSL = numpy.zeros(2*M,'d')
            fastbin.YLS = numpy.zeros(2*M,'d'); fastbin.YSL = numpy.zeros(2*M,'d')
            fastbin.ZLS = numpy.zeros(2*M,'d'); fastbin.ZSL = numpy.zeros(2*M,'d')
            
            FastGalacticBinary.FastBinaryCache[(N,T,dt)] = fastbin
        else:
            fastbin = FastGalacticBinary.FastBinaryCache[(N,T,dt)]
        
        if vector != None:
            fastbin.Response(vector[0],vector[1],0.5*math.pi - vector[2],vector[3],
                             vector[4],vector[5],vector[6],vector[7],
                             fastbin.XLS,fastbin.XSL,fastbin.YLS,fastbin.YSL,fastbin.ZLS,fastbin.ZSL)           
        else:
            fastbin.Response(self.Frequency,self.FrequencyDerivative,0.5*math.pi - self.EclipticLatitude,self.EclipticLongitude,
                             self.Amplitude,self.Inclination,self.Polarization,self.InitialPhase,
                             fastbin.XLS,fastbin.XSL,fastbin.YLS,fastbin.YSL,fastbin.ZLS,fastbin.ZSL)
        
        if buffer == None:
            # note the first frequency bin in these arrays is self.Frequency - (M/2) df, with df = 1/T
            retX, retY, retZ = map(lambda a: FrequencyArray.FrequencyArray(a[::2] - 1j * a[1::2],   # again, empirically found...
                                                                           dtype = numpy.complex128,
                                                                           kmin = int(self.Frequency*T) - M/2,df = 1.0/T),
                                   (fastbin.XSL, fastbin.YSL, fastbin.ZSL) if simulator == 'synthlisa' else (fastbin.XLS, fastbin.YLS, fastbin.ZLS))
        
            return (retX,retY,retZ)
        else:
            for i,a in enumerate((fastbin.XSL, fastbin.YSL, fastbin.ZSL) if simulator == 'synthlisa' else (fastbin.XLS, fastbin.YLS, fastbin.ZLS)):
                buffer[i][(int(self.Frequency*T) - M/2):(int(self.Frequency*T) + M/2)] += a[::2] - 1j * a[1::2]
    
    # total observation time is hardcoded to 2^22 * 15 seconds
    def TDI(self,T=6.2914560e7,dt=15.0,simulator='synthlisa',table=None):
        X, Y, Z = self.fourier(simulator,table,T=T,dt=dt)
        
        return X.ifft(dt), Y.ifft(dt), Z.ifft(dt)
    

%}