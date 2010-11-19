import re, os, sys
from distutils.dep_util import newer, newer_group
from os.path import isfile

class flow(object):
    def __init__(self,workdir="",debug=True,dryrun=False):
        self.workdir, self.debug, self.dryrun = workdir, debug, dryrun
        
        self.indent = ''
        self.makerules = {}
        self.making, self.drymade = [], []
    
    def w(self,filename):
        return filename if self.workdir == "" else self.workdir + '/' + filename
    
    def rule(self,target,sources,command):
        target  = target % self.__dict__
        sources = [s % self.__dict__ for s in sources]
        
        if type(command) == list:
            command = '; '.join(command)
        
        command = command % self.__dict__
        
        command = re.sub('\$0',self.w(target),command)
        
        for i,source in enumerate(sources):
            command = re.sub('\$%d' % (i+1),self.w(source),command)
        
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
                if isfile(self.w(target)):
                    return
                else:
                    raise ValueError, ("flow.make: Don't know how to make %s!" % target)
            
            sources,command = self.makerules[target]
            
            if self.debug: print "flow.make: %sChecking sources for target %s." % (self.indent,target); sys.stdout.flush()
            
            self.indent += '  '
            for s in sources:
                self.make(s)
            self.indent = self.indent[:-2]
            
            if (not isfile(self.w(target)) or newer_group([self.w(s) for s in sources],self.w(target))) and not (self.dryrun and target in self.drymade):
                print "flow.make: %sMaking target %s by running %s." % (self.indent,target,command); sys.stdout.flush()
                
                if self.dryrun:
                    self.drymade.append(target)
                else:
                    res = os.system(command)
                    if res != 0: sys.exit(res)
            else:
                if self.debug: print "flow.make: %sTarget %s is up to date." % (self.indent,target); sys.stdout.flush()
            
            self.making.remove(target)
    
