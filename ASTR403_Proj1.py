# ASTR 403 Project 1
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
import math



f = open("lcparam.txt", "r")
f_list = f.readlines()

# Remove labels
f_list.pop(0)

# zcmb = redshift
# mb   = observed B-band magnitude
# dmb = error of observed B-band magnitude

line = []
zcmb = []
mb = []
dmb = []

# Grab columns of data
for i in range(len(f_list)): 
    line1 = f_list[i].split(' ')
    
    zcmb.append(float(line1[1])) #Redshift
    zcmb.sort()
    mb.append(float(line1[4])) #Apparent magnitude
    mb.sort()
    dmb.append(float(line1[5])) #Error of Apparent mag. 
    dmb.sort()

# Hubble Diagram
#DM = mb - Mb = 5*log_10(DL/Mpc)+25


DM = []
Mb = -19.3
for j in range(len(mb)): 
    DM.append(float(mb[j]) - Mb)


# Plot theoretical DM 

# First model
H0 = 70*1000
Omega_m = 0.3
Omega_L = 0.7
q = 0.5*Omega_m - Omega_L
c = 3e8
#TDM = 43.23 - 5*np.log10(H0 / (68*1000)) + 5*np.log10(z) + 1.086*(1-q)*z

# Theoretical distance modulus for 1st model
TDM1 = []
for k in range(len(zcmb)): 
    
    dL = (c/H0)*zcmb[k]*(1+(1-q)*zcmb[k]/2)
    TDM_fun = 5*math.log10(dL)+25
    TDM1.append(TDM_fun)
    
    
# Second model
H0 = 70*1000
Omega_m = 1
Omega_L = 1
q = 0.5*Omega_m - Omega_L
c = 3e8
#TDM = 43.23 - 5*np.log10(H0 / (68*1000)) + 5*np.log10(z) + 1.086*(1-q)*z

# Theoretical distance modulus for 2nd model
TDM2 = []
for k in range(len(zcmb)): 
    dL = (c/H0)*zcmb[k]*(1+(1-q)*zcmb[k]/2)
    TDM_fun = 5*math.log10(dL)+25
    TDM2.append(TDM_fun)
    
    
# Third model
H0 = 70*1000
Omega_m = 0.3
Omega_L = 0
q = 0.5*Omega_m - Omega_L
c = 3e8
#TDM = 43.23 - 5*np.log10(H0 / (68*1000)) + 5*np.log10(z) + 1.086*(1-q)*z

# Theoretical distance modulus for 3rd model
TDM3 = []
for k in range(len(zcmb)): 
    dL = (c/H0)*zcmb[k]*(1+(1-q)*zcmb[k]/2)
    TDM_fun = 5*math.log10(dL)+25
    TDM3.append(TDM_fun)
    

# Find difference
Diff1 = []; Diff2 = []; Diff3 = []
for i in range(len(DM)):
    Diff1.append(DM[i] - TDM1[i])
    Diff2.append(DM[i] - TDM2[i])
    Diff3.append(DM[i] - TDM3[i])


x = symbols('x')

# plt.xscale('log')
# plt.scatter(zcmb, DM, color='blue', zorder=1)
# plt.errorbar(zcmb, DM, yerr=dmb, fmt="o")
# plt.plot(zcmb, TDM1, color='red', zorder=5, label=r'$\Omega_M = 0.3, \Omega_\Lambda = 0.7$', linewidth=3)
# plt.plot(zcmb, TDM2, color='limegreen', zorder=5, label=r'$\Omega_M = 1.0, \Omega_\Lambda = 1.0$', linewidth=2)
# plt.plot(zcmb, TDM3, color='darkorange', zorder=5, label=r'$\Omega_M = 0.3, \Omega_\Lambda = 0.0$', linewidth=2)
# plt.legend()
# plt.show()
#
# plt.xscale('log')
# plt.scatter(zcmb, Diff1)
# plt.scatter(zcmb, Diff2)
# plt.scatter(zcmb, Diff3)
# plt.ylim((-1.2,0.5))
# plt.show()
