class data(object):
    # a little magic to allow logL functions with only one parameter...
    # this is fun, but could it be done more cleanly?
    
    def setlogL(self,func):
        if func.func_code.co_argcount == 1:
            self.__dict__['logL'] = lambda signal, state: func(signal)
        else:
            self.__dict__['logL'] = func
    
    def getlogL(self):
        return self.__dict__['logL']
    
    logL = property(fget=getlogL,fset=setlogL)

