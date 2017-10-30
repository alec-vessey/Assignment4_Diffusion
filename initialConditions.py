import numpy as np

# Initial conditions function for diffusion

def squareWave(x,alpha,beta):
    "A square wave as a function of position, x, which is 1 between alpha"
    "and beta and zero elsewhere. The initialisation is conservative so"
    "that each phi contains the correct quantity integrated over a region"
    "a distance dx/2 either side of x"
    
    phi = np.zeros_like(x)
    
    # The grid spacing (assumed uniform)
    dx = x[1] - x[0]
    #print dx
    # Set phi away from the end points (assume zero at the end points)
    for j in xrange(1,len(x)-1):    #1,len(x)-1
        # edges of the grid box (using west and east notation)
        xw = x[j] - 0.5*dx
        xe = x[j] + 0.5*dx
        
        #integral quantity of phi
        phi[j] = max((min(beta, xe) - max(alpha, xw))/dx, 0)
        
        
        #print (dx)
        #print (min(beta, xe), max(alpha, xw), dx)
        
        #print max(0.0, 0.4)

    return phi

def squareWave_alternative(x,alpha,beta):
    "A square wave as a function of position, x, which is 1 between alpha"
    "and beta and zero elsewhere. The initialisation is conservative so"
    "that each phi contains the correct quantity integrated over a region"
    "a distance dx/2 either side of x"
    
    phi = np.zeros_like(x)
    
    # The grid spacing (assumed uniform)
    dx = x[1] - x[0]
    #print dx
    # Set phi away from the end points (assume zero at the end points)
    for j in xrange(1,len(x)-1):    #1,len(x)-1
        # edges of the grid box (using west and east notation)
#        xw = x[j] - 0.5*dx
#        xe = x[j] + 0.5*dx
#        
#        #integral quantity of phi
#        phi[j] = max((min(beta, xe) - max(alpha, xw))/dx, 0)
#        
#        #print (dx)
#        #print (min(beta, xe), max(alpha, xw), dx)
#        
#        #print max(0.0, 0.4)
        #print (j*dx)
        if j*dx >= 0.4 and j*dx < 0.625:
            phi[j] = 1
        else:
            phi[j] = 0

    return phi

xmin = 0 
xmax = 1 
nx = 41
squareWaveMin = 0.4
squareWaveMax = 0.6

dx = (xmax - xmin)/(nx-1)

x = np.zeros(nx)
for j in xrange(nx):
    x[j] = xmin + j*dx

#print x
#squareWaveMin = 0.4
#squareWaveMax = 0.6
#
#print squareWave(x, squareWaveMin, squareWaveMax)