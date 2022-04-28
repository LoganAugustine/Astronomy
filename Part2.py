# ASTR 403 Project 1
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
import math
import emcee

import warnings
warnings.filterwarnings('ignore')


# Ready? Begin.
f = open("lcparam.txt", "r")
f_list = f.readlines()

# Remove labels
f_list.pop(0)

# zcmb = redshift
# mb   = observed B-band magnitude
# dmb = error of observed B-band magnitude

# Grab columns of data and append
line = []
zcmb = []
mb = []
dmb = []
for i in range(len(f_list)):
    line1 = f_list[i].split(' ')

    zcmb.append(float(line1[1]))  # Redshift
    zcmb.sort()
    mb.append(float(line1[4]))  # Apparent magnitude
    mb.sort()
    dmb.append(float(line1[5]))  # Error of Apparent mag.
    dmb.sort()


# Section 1 on Hubble Diagram

DM = []
Mb = -19.3
for j in range(len(mb)):
    DM.append(float(mb[j]) - Mb)


# Plot theoretical DM

# First model
H0 = 70 * 1000
Omega_m = 0.3
Omega_L = 0.7
q = 0.5 * Omega_m - Omega_L
c = 3e8
# DM formula = mb - Mb = 5*log_10(DL/Mpc)+25
# TDM formula = 43.23 - 5*np.log10(H0 / (68*1000)) + 5*np.log10(z) + 1.086*(1-q)*z

# TDM for 1st model
TDM1 = []
for k in range(len(zcmb)):
    dL = (c / H0) * zcmb[k] * (1 + (1 - q) * zcmb[k] / 2)
    TDM_fun = 5 * math.log10(dL) + 25
    TDM1.append(TDM_fun)

# Second model
H0 = 70 * 1000
Omega_m = 1
Omega_L = 1
q = 0.5 * Omega_m - Omega_L
c = 3e8
# TDM = 43.23 - 5*np.log10(H0 / (68*1000)) + 5*np.log10(z) + 1.086*(1-q)*z

# TDM for 2nd model
TDM2 = []
for k in range(len(zcmb)):
    dL = (c / H0) * zcmb[k] * (1 + (1 - q) * zcmb[k] / 2)
    TDM_fun = 5 * math.log10(dL) + 25
    TDM2.append(TDM_fun)

# Third model
H0 = 70 * 1000
Omega_m = 0.3
Omega_L = 0
q = 0.5 * Omega_m - Omega_L
c = 3e8
# TDM = 43.23 - 5*np.log10(H0 / (68*1000)) + 5*np.log10(z) + 1.086*(1-q)*z

# TDM for 3rd model
TDM3 = []
for k in range(len(zcmb)):
    dL = (c / H0) * zcmb[k] * (1 + (1 - q) * zcmb[k] / 2)
    TDM_fun = 5 * math.log10(dL) + 25
    TDM3.append(TDM_fun)

# Find difference
Diff1 = []
Diff2 = []
Diff3 = []
for i in range(len(DM)):
    Diff1.append(DM[i] - TDM1[i])
    Diff2.append(DM[i] - TDM2[i])
    Diff3.append(DM[i] - TDM3[i])

x = symbols('x')

DM_A = np.array(DM)
zcmb_A = np.array(zcmb)
dmb_A = np.array(dmb)
yerr = dmb_A
x = zcmb_A
y = DM_A
H0 = 70
OM = 0.3
x0 = np.linspace(0.01, 2.27, 500)

from astropy.cosmology import FlatLambdaCDM

def DMFORM(x, H0, OM):
    OL = 1 - OM
    q0 = 0.5 * OM - OL
    c = 3e8
    dl = (c/H0)*x*(1 + (1-q0)*x/2)

    cosmo = FlatLambdaCDM(H0=H0, Om0=OM)

    Dl = cosmo.luminosity_distance(x).value

    # finding distance modulus
    dist_mod = 5 * np.log10(Dl) + 25

    return (dist_mod)


# Begin emcee tutorial
def log_likelihood(theta, x, y, yerr):
    H0, OM = theta
    model = DMFORM(x, H0, OM)
    sigma2 = yerr ** 2
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))


from scipy.optimize import minimize

np.random.seed(42)
nll = lambda *args: -log_likelihood(*args)
initial = np.array([H0, OM])
soln = minimize(nll, initial, args=(x, y, yerr))
H0_ml, OM_ml = soln.x

print("Maximum likelihood estimates:")
print("H0 = {0:.3f}".format(H0_ml))
print("OM = {0:.3f}".format(OM_ml))
print("OL = {0:.3f}".format(1-OM_ml))


plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)
plt.plot(x, DMFORM(x,H0_ml,OM_ml), "k", alpha=0.3, lw=3, label="truth")
# plt.plot(x0, DMFORM(x0, fit[0], fit[1], "--k", label="LS")
plt.legend(fontsize=14)
plt.xlabel("x")
plt.ylabel("y")
plt.show()


def log_prior(theta):
    H0, OM = theta
    if 60. < H0 < 80. and 0. < OM < 1.:
        return 0.0
    return -np.inf


def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)


pos = soln.x + 1e-4 * np.random.randn(32,2)

nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, log_probability, args=(x, y, yerr)
)
sampler.run_mcmc(pos, 5000, progress=True);


fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["m", "b", "log(f)"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number")
plt.show()


tau = sampler.get_autocorr_time()
print(tau)


flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
print(flat_samples.shape)


import corner

fig = corner.corner(
    flat_samples, labels=labels, truths=[H0, OM]
)
plt.show()


inds = np.random.randint(len(flat_samples), size=100)
for ind in inds:
    sample = flat_samples[ind]
    plt.plot(x0, DMFORM(x0, sample[0], sample[1]), "C1", alpha=0.1)
plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)

plt.plot(x0, DMFORM(x0, sample[0], sample[1]), "k", label="truth")
plt.legend(fontsize=14)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

emcee_vals = []

for i in range(ndim):
    mcmc = np.percentile(flat_samples[i:, i], [16, 50, 84])
    emcee_vals.append(mcmc[1])



print('H0:', emcee_vals[0])
print('OM:', emcee_vals[1])

# H0 = 71, OM = 0.28



