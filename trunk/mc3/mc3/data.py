from utility import stripargument

class data(object):
    def setlogL(self,func):
        self.__dict__['logL'] = stripargument(func)
    
    def getlogL(self):
        return self.__dict__['logL']
    
    logL = property(fget=getlogL,fset=setlogL)

