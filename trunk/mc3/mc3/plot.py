import math

import matplotlib.pyplot as P
import numpy as N

def normaldist(mean,var):
    return lambda x: N.exp(-0.5*(x-mean)**2/var) / math.sqrt(2.0*math.pi*var)

def plotchain(chain,title='MCMC parameter plots',fignum=1,otherpars=[],xbins=50,corrwindow=100,overplot=None,plotbins=100,index=0):
    P.figure(fignum)
    
    plotpars = [par for par in chain.model.parameters if par in chain.tostore]
    
    subplots = len(plotpars) + len(otherpars)
        
    # replicate plot parameter values for all parameters if not given explicitly
    xbins, corrwindow, overplot, plotbins = (par if isinstance(par,(tuple,list)) else [par] * subplots
                                             for par in (xbins, corrwindow, overplot, plotbins))
    
    for i,par in enumerate(plotpars + otherpars):
        x = N.array(chain[par])
        
        if x.ndim > 1:  # if chain consists of multistates, choose a section
            x = x[:,index]
        
        # make a histogram
        P.subplot(subplots,3,3*i+1); P.hist(x,bins=xbins[i],normed=True); P.xlabel(par)
        if overplot[i]:
            a = N.linspace(N.min(x),N.max(x),plotbins[i])
            P.plot(a,overplot[i](a))
        # make a trajectory plot
        P.subplot(subplots,3,3*i+2); P.plot(x); P.xlabel(par)
        # make a normalized correlation plot
        xm = x - N.mean(x)
        c = N.correlate(xm,xm,'full')[len(x)-1:len(x)+corrwindow[i]]; c = c / c[0]
        P.subplot(subplots,3,3*i+3); P.plot(c); P.xlabel(par)
    
    P.suptitle(title)

def statchain(chain,report=True,skip=1,cov=True,err=False):
    if report:
        total = len(chain.samples)
        accepted = sum(1 for i in range(1,total) if chain.samples[i] != chain.samples[i-1])
        
        print "Chain report: %d samples, %d moves accepted (%.2f)" % (total,accepted,1.0*accepted/total)
        
    parray = N.array(chain.samples[::skip],'d')
    if skip > 1:
        print "(Using every %dth sample:)" % skip
        
    mean = N.mean(parray,axis=0)
    print "Cond. means :" + (" %g" * len(mean)) % tuple(mean)
    
    covmat = N.cov(parray.T)
    
    if cov:
        print "Covariance  :"; print covmat
    
    if err:
        print "Errors      :"; print N.sqrt(N.diagonal(covmat))
    
