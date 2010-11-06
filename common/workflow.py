import re, os, sys
from distutils.dep_util import newer, newer_group
from os.path import isfile

class flow(object):
    def __init__(self,workdir="",commands=None,debug=True,dryrun=False):
        self.debug, self.dryrun = debug, dryrun
        if self.dryrun: self.drymade = []
        
        if workdir != "": self.workdir = workdir + '/'        
        
        self.makerules = {}
        self.indent = ''
        self.making = []
    
    def rule(self,target,sources,command):
        target  = target % self.__dict__
        sources = [s % self.__dict__ for s in sources]
        command = command % self.__dict__
        
        command = re.sub('\$0',self.workdir + target,command)
        
        for i,source in enumerate(sources):
            command = re.sub('\$%d' % (i+1),self.workdir + source,command)
        
        self.makerules[target] = (sources,command)
    
    def make(self,target=None):
        if target == None:
            for t in self.makerules:
                self.make(t)            # these won't be in any particular order, but it does not matter
        else:
            if target in self.making:
                raise RuntimeError, ("flow.make: Met a circular condition while making %s!" % target)
            else:
                self.making.append(target)
            
            if target not in self.makerules:
                if isfile(self.workdir + target):
                    return
                else:
                    print self.makerules
                    raise ValueError, ("flow.make: Don't know how to make %s!" % target)
            
            sources,command = self.makerules[target]
            
            if self.debug: print "flow.make: %sChecking sources for target %s." % (self.indent,target); self.indent += '  '
            for s in sources:
                self.make(s)
            if self.debug: self.indent = self.indent[:-2]
                    
            if (not isfile(self.workdir + target) or newer_group([self.workdir + s for s in sources],self.workdir + target)) and not (self.dryrun and target in self.drymade):
                if self.debug: print "flow.make: %sMaking target %s." % (self.indent,target)
                
                if self.dryrun:
                    print "flow.make: RUNNING %s" % command
                    self.drymade.append(target)
                else:
                    # print "flow.make: RUNNING %s" % command
                    # os.system('touch %s%s' % (self.workdir,target)) 
                    res = os.system(command)
                    if res != 0: sys.exit(res)
            else:
                if self.debug: print "%sTarget %s is up to date." % (self.indent,target)
            
            self.making.remove(target)
    
