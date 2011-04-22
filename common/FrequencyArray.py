import numpy

class FrequencyArray(numpy.ndarray):
    def __new__(subtype,data,dtype=None,copy=False,df=None,kmin=None):
        # make sure we are working with an array, copy the data if requested,
        # then transform the array to our new subclass
        subarr = numpy.array(data,dtype=dtype,copy=copy)
        subarr = subarr.view(subtype)
        
        # get df and kmin preferentially from the initialization,
        # then from the data object, otherwise set to None
        subarr.df   = df   if df   is not None else getattr(data,'df',  None)
        subarr.kmin = kmin if kmin is not None else getattr(data,'kmin',0)
        
        return subarr
    
    def __array_wrap__(self,out_arr,context=None):
        out_arr.df, out_arr.kmin = self.df, self.kmin
        
        return numpy.ndarray.__array_wrap__(self,out_arr,context)
        
    def __getitem__(self,key):
        return self.view(numpy.ndarray)[key]
    
    def __getslice__(self,i,j):
        return self.view(numpy.ndarray)[i:j]
    
    # def __array_finalize__(self,obj):
    #    if obj is None: return
    #    self.df   = getattr(obj,'df',  None)
    #    self.kmin = getattr(obj,'kmin',None)
    
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
            
            ret = numpy.zeros(end-beg,dtype=self.dtype)
            
            ret[(self.kmin  - beg):(self.kmin  - beg + len(self))]   = self
            ret[(other.kmin - beg):(other.kmin - beg + len(other))] += other
            
            return FrequencyArray(ret,kmin=beg,df=self.df)
        
        # fall back to simple arrays (may waste memory)
        return numpy.ndarray.__add__(self,other)
    
    # same behavior as __add__: TO DO -- consider restricting the result to the extend of the first array
    def __sub__(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            beg = min(self.kmin,other.kmin)
            end = max(self.kmin + len(self),other.kmin + len(other))
            
            ret = numpy.zeros(end-beg,dtype=self.dtype)
            
            ret[(self.kmin  - beg):(self.kmin  - beg + len(self))]   = self
            ret[(other.kmin - beg):(other.kmin - beg + len(other))] -= other
            
            return FrequencyArray(ret,kmin=beg,df=self.df)
        
        # fall back to simple arrays (may waste memory)
        return numpy.ndarray.__sub__(self,other)        
    
    # restrict the result to the extent of the first array (useful, e.g., for logL over frequency-limited data)
    def rsub(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            if other.kmin >= self.kmin + len(self) or self.kmin >= other.kmin + len(other):
                return self
            else:
                beg = max(self.kmin,other.kmin)
                end = min(self.kmin + len(self),other.kmin + len(other))
                
                ret = numpy.array(self,copy=True)
                ret[(beg - self.kmin):(end - self.kmin)] -= other[(beg - other.kmin):(end - other.kmin)]
                
                return FrequencyArray(ret,kmin=self.kmin,df=self.df)
        
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
                ret = numpy.array(self[(beg - self.kmin):(end - self.kmin)],copy=True)
                ret *= other[(beg - other.kmin):(end - other.kmin)]
                
                return FrequencyArray(ret,kmin=beg,df=self.df)
        
        # fall back to simple arrays (may waste memory)
        return numpy.ndarray.__mul__(self,other)
    
    # in division, it's OK if second array is larger, but not if its smaller (which implies division by zero!)
    def __div__(self,other):
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            if (other.kmin > self.kmin) or (other.kmin + len(other) < self.kmin + len(self)):
                raise ZeroDivisionError
            else:
                ret = numpy.array(self,copy=True)
                ret /= other[(self.kmin - other.kmin):(self.kmin - other.kmin + len(self))]
            
            return FrequencyArray(ret,kmin=self.kmin,df=self.df)
        
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
