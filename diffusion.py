#!/usr/bin/python

# Outer code for setting up the diffusion problem on a uniform
# grid and calling the function to perform the diffusion and plot.

from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
execfile("diffusionSchemes.py")
execfile("diagnostics.py")
execfile("initialConditions.py")


def main():
    """
    Diffuse a squareWave between squareWaveMin and squareWaveMax on a
    domian between x = xmin and x = xmax split over nx spatial steps
    with diffusion coefficient K, time step dt for nt time steps
    """
    # Parameters
    xmin = 0 
    xmax = 1 
    nx = 41
    nt = 40
    dt = 0.1
    k = 1e-3
    squareWaveMin = 0.4
    squareWaveMax = 0.6
    
    # Derived parameters
    dx = (xmax - xmin)/(nx-1)
    d = k*dt/dx**2     #Non-dimensional diffusion coefficent (d)
    print ("non-dimensional diffusion coefficient = ", d)
    print ("dx = ", dx, "dt = ", dt, "nt = ", nt)
    print ("end time = ", nt*dt)
    
    # spatial points for plotting and for defining initial conditions
    x = np.zeros(nx)
    for j in xrange(nx):
        x[j] = xmin + j*dx

    
    # initial conditions
    phiOld = squareWave(x, squareWaveMin, squareWaveMax)
    # analytical solution (of square wave profile in an infinite domain)
    phiAnalytic = analyticErf(x, k*dt*nt, squareWaveMin, squareWaveMax)
    
    #Diffusion using FTCS and BTCS
    phiFTCS = FTCS(phiOld.copy(), d, nt)
    phiBTCS = BTCS(phiOld.copy(), d, nt)
    
    #calculate and print out error norms
    #print("FTCS L2 error norm =", L2ErrorNorm(phiFTCS, phiAnalytic))
    #print("BTCS L2 error norm =", L2ErrorNorm(phiBTCS, phiAnalytic))
#    
#
    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color = 'black')
    plt.plot(x, phiAnalytic, label='Analytical', color='black', linestyle='--', linewidth=2)
    plt.plot(x, phiFTCS, label='FTCS', color='blue')
    plt.plot(x, phiBTCS, label='BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([0,1])
    plt.xlabel('$x$')
    plt.ylabel('$\phi$')
    plt.title('Graph showing the differences between a forward-time centred \n space (FTCS) and a backward-time centred space (BTCS) \n finite differences scheme of diffusion after 40 timesteps', y = 1.1)
    plt.legend(bbox_to_anchor=(1.4, 0.61))    
    #plt.savefig('/msc_modules/Intro_Numerical_Modelling_git/Assignment4_Diffusion/figures/figure_1a_FTCSBTCS_comp.pdf', format = 'pdf', dpi = 1000, bbox_inches='tight')    #plot the solutions
    
    #calculate the errors
    phiFTCS_error = []
    phiBTCS_error = []
    for q in xrange(len(phiFTCS)):
        phiFTCS_error.append(phiFTCS[q] - phiAnalytic[q])
        phiBTCS_error.append(phiBTCS[q] - phiAnalytic[q])
        
    #find the absolute error
    phiFTCS_error = np.abs(phiFTCS_error)
    phiBTCS_error = np.abs(phiBTCS_error)
    
    #plot the errors
    plt.figure(2)
    plt.plot(x, phiFTCS_error, label = 'FTCS Absolute Error', color = 'blue')
    plt.plot(x, phiBTCS_error, label = 'BTCS Absolute Error', color = 'red')
    plt.title('Graph showing the absolute error of a forward-time centred space (FTCS) \n and a backward-time centred space (BTCS) finite differences \n scheme of diffusion', y = 1.1)    
    plt.ylim(ymax = 0.05)    
    plt.xlabel('$x$')
    plt.ylabel('$\phi$ Error')
    plt.legend(bbox_to_anchor=(1.55, 0.6))
    #plt.show()
    ##plt.savefig('/msc_modules/Intro_Numerical_Modelling_git/Assignment4_Diffusion/figures/figure_1b_errors.pdf', format = 'pdf', dpi = 1000, bbox_inches='tight')
#    
##    """Question 3"""
    # run analystical and FTCS for long duration so that it gives significantly different results
    nt = 1000
    phiAnalytic_long_duration = analyticErf(x, k*dt*nt, squareWaveMin, squareWaveMax)
    phiFTCS_long_duration = FTCS(phiOld.copy(), d, nt)
    
    plt.figure(3)
    plt.plot(x, phiAnalytic_long_duration, label = 'Analytical Long \n Duration', color = 'black', linestyle='--')
    plt.plot(x, phiFTCS_long_duration, label = 'FTCS Long \n Duration', color = 'blue')
    plt.title('Graph showing the differences between a forward time centred space (FTCS) \n finite differences scheme and the analystical solution over \n a long duration (1000 time-steps)', y = 1.1)       
    plt.xlabel('$x$')
    plt.ylabel('$\phi$')
    plt.xlim(xmin = 0.0)
    plt.xlim(xmax = 1.0)
    plt.legend(bbox_to_anchor=(1.55, 0.6))
    #plt.show()
    plt.savefig('/msc_modules/Intro_Numerical_Modelling_git/Assignment4_Diffusion/figures/figure_2_longn_comp.pdf', format = 'pdf', dpi = 1000, bbox_inches='tight')
#    
#    """Question 4"""
#    #FTCS stability
#    #for FTCS d (the non-dimensional diffusion coefficient) needs to be less than 0.5
#    #you expect the error to increase when d becomes higher than 0.5 - plot grpahs for 
#    #different values of d and see if the value increases     
#    #change d by changing the values of dt
#    
#    # set parameters
#    ##increase the time steps
    nt = [40, 80, 120, 160] 
    graph_colors = ['blue', 'red', 'green', 'black']
#    
#    ##change the value of k so that d is greater than 0.5, hence making FTCS 
#    ##unstable
    k = 5e-3
    d = k*dt/dx**2 
    
    #plot the graph of errors depending on the value of d
    plt.figure(4)
    d_array = []

    for i in xrange(len(nt)):
        #calculate d then calculate FTCS
        phiFTCS = FTCS(phiOld.copy(), d, nt[i])
        phiAnalytic = analyticErf(x, k*dt*nt[i], squareWaveMin, squareWaveMax)
        #calculate the error
        phiFTCS_error = []
        L2ErrorNorm_error = []
        for q in xrange(len(phiFTCS)):
            phiFTCS_error.append(phiFTCS[q] - phiAnalytic[q])
        L2ErrorNorm_error.append(L2ErrorNorm(phiFTCS, phiAnalytic))
        #find the absolue error
        phiFTCS_error = np.abs(phiFTCS_error)
        print ('l2error norm tendency:', L2ErrorNorm_error)
        #plot the error
        plt.plot(x, phiFTCS_error, label = 'nt = {}'.format(nt[i]), color = graph_colors[i])
    
    plt.title('Graph showing the absolute error of a forward-time centred space (FTCS) finite differences \n scheme with varying lengths of time when d is equal to {}'.format(d), y = 1.1)      
    plt.xlabel('$x$')
    plt.ylabel('$\phi$ Error')
    plt.legend(bbox_to_anchor=(1.55, 0.6))
    #plt.show()
    plt.savefig('/msc_modules/Intro_Numerical_Modelling_git/Assignment4_Diffusion/figures/figure_3a_errors_FTCS_d016.pdf', format = 'pdf', dpi = 1000, bbox_inches='tight')    
#    
#    #BTCS stability 
#    #for FTCS d (the non-dimensional diffusion coefficient) needs to be less than 0.0
#    #but d can never be negative, therefore, BTCS is unconditionally stable   

    
    """Question 5"""
    # a finite difference scheme is convergent if the solution of the scheme 
    ##converge to the solution of the PDE as dx and dt converge to 0
#    d = k*dt/dx**2 
    #spatial resolution
    xmin = 0
    xmax = 1
    nx_string = [21, 41, 81]
    k = 1e-3
    d = 0.16
    squareWaveMin = 0.4
    squareWaveMax = 0.6
    l2FTCS_string = []
    dx_string = []
    error_string = []
    dt_string = []
    
    #calculate phiold, phiFTCS and phiAnalytical for each change in dt and nx
    for n in range(0, 3, 1):
        dx = (xmax - xmin)/(nx_string[n]-1)
        dx_string.append(dx)
        dt = d*(dx**2)/k 
        nt = 4.0/dt
        dt_string.append(dt)
        print ('|dx|dt|nt|nt*dt|d|')
        print (dx, dt, nt, nt*dt, d)
        x = np.zeros(nx_string[n])
        
        for j in xrange(nx_string[n]):
            x[j] = xmin + j*dx
    
        phiOld = squareWave(x, squareWaveMin, squareWaveMax)
        phiFTCS = FTCS(phiOld.copy(), d, nt)
        phiAnalytic = analyticErf(x, k*dt*nt, squareWaveMin, squareWaveMax)
        
        #calculate l2 error and append to array
        l2FTCS_string.append(L2ErrorNorm(phiFTCS, phiAnalytic))
        
        #calculate absolute error
        error = np.abs(phiFTCS[0] - phiAnalytic[0])
        error_string.append(error)
        
    print ('error:', error_string)
    
    # plot l2 error and dx
    plt.figure(6)
    plt.plot(dx_string, l2FTCS_string, label = 'error')
    plt.scatter(dx_string, l2FTCS_string, label = 'error')
    plt.xlabel('dt')
    plt.ylabel('L2FTCS')
    plt.xticks(np.arange(0.02, 0.4, 0.06))
    #plt.yticks(np.arange(min(l2FTCS_string), max(l2FTCS_string), 1))
    #plt.yscale('log')
    #plt.xscale('log')
    plt.legend()
    plt.show()
    #plt.savefig('/Users/Alec/Desktop/l2_error_dt.pdf', format = 'pdf', dpi = 1000, bbox_inches='tight')
    
#    slope, intercept = np.polyfit(np.log(dt_string), np.log(error_string), 1)
#    print('slope:', slope)

    #calculate the order of convergence
    FTCS_order_of_accuracy = (np.log(error_string[0]) - np.log(error_string[1])) / (np.log(dt_string[0]) - np.log(dt_string[1]))
    print ('order of convergence:', FTCS_order_of_accuracy)

    """ Question 6 """
#    
    xmin = 0 
    xmax = 1 
    nx = 41
    nt = 40
    dt = 0.1
    k = 1e-3
    squareWaveMin = 0.4
    squareWaveMax = 0.6
    
    # Derived parameters
    dx = (xmax - xmin)/(nx-1)
    d = k*dt/dx**2     #Non-dimensional diffusion coefficent (d)
    print ("non-dimensional diffusion coefficient = ", d)
    print ("dx = ", dx, "dt = ", dt, "nt = ", nt)
    print ("end time = ", nt*dt)
    
    # spatial points for plotting and for defining initial conditions
    x = np.zeros(nx)
    for j in xrange(nx):
        x[j] = xmin + j*dx
       
    phiOld = squareWave(x, squareWaveMin, squareWaveMax)
    phiOld_alternative = squareWave_alternative(x, squareWaveMin, squareWaveMax)
    
    phiAnalytic = analyticErf(x, k*dt*nt, squareWaveMin, squareWaveMax)
    
    #DIffusion using FTCS and BTCS
    phiFTCS = FTCS(phiOld.copy(), d, nt)
    phiBTCS = BTCS(phiOld.copy(), d, nt)
    
    #print (phiOld, phiOld_alternative)
    phiFTCS_alternative = FTCS(phiOld_alternative.copy(), d, nt)    
    phiBTCS_alternative =BTCS(phiOld_alternative.copy(), d, nt)
    
    #calculate error of initial scheme
    L2_FTCS = L2ErrorNorm(phiFTCS, phiAnalytic)
    L2_BTCS = L2ErrorNorm(phiBTCS, phiAnalytic)
    
    #calcuate error of alternative scheme
    L2_FTCS_alternative = L2ErrorNorm(phiFTCS_alternative, phiAnalytic)
    L2_BTCS_alternative =  L2ErrorNorm(phiBTCS_alternative, phiAnalytic)
    
    #print (L2_FTCS, L2_FTCS_alternative)
    #print (L2_BTCS, L2_BTCS_alternative)
    
    #print (L2_FTCS, L2_FTCS_alternative)
    #print (L2_BTCS, L2_BTCS_alternative)
    
    font = {'size' : 10}    
    plt.rc('font', **font)
    plt.figure(7)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial Conditions', color = 'black')
    plt.plot(x, phiOld_alternative, label='Alternative \n Initial Conditions', color = 'black', linestyle='--')
    #plt.plot(x, phiAnalytic, label='Analytical', color='black', linestyle='--', linewidth=2)
    #plt.plot(x, phiFTCS, label='FTCS', color='blue')
    #plt.plot(x, phiFTCS_alternative, label='FTCS', color='blue', linestyle='--')
    #plt.plot(x, phiBTCS, label='BTCS', color='red')
    #plt.plot(x, phiBTCS_alternative, label='BTCS', color='red', linestyle='--')
    plt.ylim([0,1])
    plt.xlim([0,1])
    plt.xlabel('$x$')
    plt.ylabel('$\phi$')
    #plt.title('Graph showing the differences between the initial conditions when the jump from 0.0 to 1.0 is', y = 1.1)
    plt.legend(bbox_to_anchor=(1.5, 0.61))    
    #plt.show()
    plt.savefig('/msc_modules/Intro_Numerical_Modelling_git/Assignment4_Diffusion/figures/figure_5_initial_conditions.pdf', format = 'pdf', dpi = 1000, bbox_inches='tight')

main()
    
    
