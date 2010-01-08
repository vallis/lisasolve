import math

class state(dict):
    def __init__(self,model):
        self.model = model

    # it may be good to cache prior
    
    def getprior(self):
        return self.model.prior(self)
    
    prior = property(getprior)
    
    # the idea here is to use access parameters defined in model.extra    
    def __getitem__(self,attr):
        if attr in self:
            return dict.__getitem__(self,attr)
        else:
            return self.model.extra[attr]
    
