import numpy

class FrequencyArray(numpy.ndarray):
    defaults = {'df': None, 'kmin': 0}
    
    def __new__(subtype, data, dtype=None, copy=False, **kwargs):
        # Make sure we are working with an array, and copy the data if requested
        subarr = numpy.array(data, dtype=dtype, copy=copy)
        
        # Transform 'subarr' from an ndarray to our new subclass.
        subarr = subarr.view(subtype)
        
        for attr in ['df','kmin']:
            if (attr in kwargs):
                setattr(subarr,attr,kwargs[attr])
            elif hasattr(data,attr):
                setattr(subarr,attr,getattr(data,attr))
        
        # Finally, we must return the newly created object:
        return subarr
    
    def __array_finalize__(self,obj):
        # this sets the defaults
        for attr in ['df','kmin']:      
            if not hasattr(self, attr):
                setattr(self,attr,getattr(obj,attr,self.defaults[attr]))
    
    def __repr__(self):
        if self.df is not None:
            return 'Frequency array (f0=%s,df=%s): %s' % (self.kmin * self.df,self.df,self)
        else:
            return str(self)
    
    # combine two FrequencyArrays into a longer one by adding intermediate zeros if necessary
    def __add__(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            beg = min(self.kmin,other.kmin)
            end = max(self.kmin + len(self),other.kmin + len(other))
            
            ret = FrequencyArray(numpy.zeros(end-beg,dtype=self.dtype),kmin=beg,df=self.df)
            
            ret[(self.kmin  - ret.kmin):(self.kmin  - ret.kmin + len(self) )]  = self[:]
            numpy.ndarray.__iadd__(ret[(other.kmin - ret.kmin):(other.kmin - ret.kmin + len(other))],other[:])
            
            return ret
        
        # fall back to simple arrays (may waste memory)
        return numpy.ndarray.__add__(self,other)
    
    def __sub__(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            beg = min(self.kmin,other.kmin)
            end = max(self.kmin + len(self),other.kmin + len(other))
            
            ret = FrequencyArray(numpy.zeros(end-beg,dtype=self.dtype),kmin=beg,df=self.df)
            
            ret[(self.kmin  - ret.kmin):(self.kmin  - ret.kmin + len(self) )]  = self[:]
            numpy.ndarray.__isub__(ret[(other.kmin - ret.kmin):(other.kmin - ret.kmin + len(other))],other[:])
            
            return ret
        
        # fall back to simple arrays (may waste memory)
        return numpy.ndarray.__sub__(self,other)        
    
    # rsub will restrict the result to the extent of the first array,
    # which makes sense for data, if you think about it
    def rsub(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            if other.kmin >= self.kmin + len(self) or self.kmin >= other.kmin + len(other):
                return self
            else:
                beg = max(self.kmin,other.kmin)
                end = min(self.kmin + len(self),other.kmin + len(other))
                
                ret = self.copy()
                
                numpy.ndarray.__isub__(ret[(beg - self.kmin):(end - self.kmin)],other[(beg - other.kmin):(end - other.kmin)])
                
                return ret
        
        return numpy.ndarray.__sub__(self,other)
    
    # the inplace add and sub will work only if the second array is contained in the first one
    def __iadd__(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            if (self.kmin <= other.kmin) and (self.kmin + len(self) >= other.kmin + len(other)):
                numpy.ndarray.__iadd__(self[(other.kmin - self.kmin):(other.kmin - self.kmin + len(other))],other[:])
                return self
        
        # fall back to simple arrays
        numpy.ndarray.__iadd__(self,other)
        return self
    
    def __isub__(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            if (self.kmin <= other.kmin) and (self.kmin + len(self) >= other.kmin + len(other)):
                numpy.ndarray.__isub__(self[(other.kmin - self.kmin):(other.kmin - self.kmin + len(other))],other[:])
                return self
        
        # fall back to simple arrays
        numpy.ndarray.__isub__(self,other)
        return self
    
    
    # in multiplication, we go for the intersection of arrays (not their union!)
    # no intersection return a scalar 0
    def __mul__(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            beg = max(self.kmin,other.kmin)
            end = min(self.kmin + len(self),other.kmin + len(other))
            
            if beg >= end:
                return 0.0
            else:
                ret = FrequencyArray(numpy.zeros(end-beg,dtype=self.dtype),kmin=beg,df=self.df)
            
                ret[:]  = self[(ret.kmin - self.kmin):(ret.kmin - self.kmin + len(ret))]
                numpy.ndarray.__imul__(ret[:],other[(ret.kmin - other.kmin):(ret.kmin - other.kmin + len(ret))])
            
                return ret
        
        # fall back to simple arrays (may waste memory)
        return numpy.ndarray.__mul__(self,other)
    
    # the inplace multiply works only if the second array is contained in the first
    def __imul__(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            if (self.kmin <= other.kmin) and (self.kmin + len(self) >= other.kmin + len(other)):
                numpy.ndarray.__imul__(self[(other.kmin - self.kmin):(other.kmin - self.kmin + len(other))],other[:])
                return self
                
        # fall back to simple arrays
        numpy.ndarray.__imul__(self,other)
        return self
    
    # division works only if the first array is contained in the second
    def __div__(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            if (self.kmin >= other.kmin) and (self.kmin + len(self) <= other.kmin + len(other)):
                ret = self.copy()
                
                numpy.ndarray.__idiv__(ret[:],other[(self.kmin - other.kmin):(self.kmin - other.kmin + len(self))])
                
                return ret
        
        # fall back to simple arrays
        return numpy.ndarray.__div__(self,other)
    
    @property
    def f(self):
        return numpy.linspace(self.kmin * self.df,(self.kmin + len(self) - 1) * self.df,len(self))
    
    @property
    def fmin(self):
        return self.kmin * self.df
    
    @property
    def fmax(self):
        return (self.kmin + len(self)) * self.df
    
    def ifft(self,dt):
        n = int(1.0/(dt*self.df))
        
        ret = numpy.zeros(n/2+1,dtype=self.dtype)
        ret[self.kmin:self.kmin+len(self)] = self[:]
        ret *= n                                        # normalization, ehm, found empirically
        
        return numpy.fft.irfft(ret)
    
    # restrict the array to the dimensions of the second, or to dimensions specified as (kmin,len)
    def restrict(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            kmin, length = other.kmin, len(other)
        elif isinstance(other,(list,tuple)) and len(other) == 2:
            kmin, length = other
        else:
            raise TypeError
        
        # no need to restrict anything?
        if kmin == self.kmin and length == len(self):
            return other
        
        ret = FrequencyArray(numpy.zeros(length,dtype=self.dtype),kmin=kmin,df=self.df)
        
        beg = max(self.kmin,kmin)
        end = min(self.kmin + len(self),kmin + length)
        
        ret[(beg - kmin):(end - kmin)] = self[(beg - self.kmin):(end - self.kmin)]
        
        return ret
    
    # pad the array on both sides
    def pad(self,leftpad=1,rightpad=1):
        return self.restrict((self.kmin - int(leftpad)*len(self),int(1+leftpad+rightpad)*len(self)))
    


# a = FrequencyArray([1,2,3,4,5],kmin=1)
# b = FrequencyArray([1,2,3,4],kmin=2)

# print 2 * a, type(a), (2*a).kmin
