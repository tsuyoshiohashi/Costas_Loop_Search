#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
# costas_loop_search.py
#
# Reference : https://pysdr.org/content/rds.html
#
import numpy as np
from scipy.signal import resample_poly, firwin
import matplotlib.pyplot as plt
import sys

# Read in IQ signal from file
#x = np.fromfile('./433M_bpsk64k_1024kSample_ran.iq', dtype=np.complex64)
args = sys.argv
if(len(sys.argv) == 2):
    path = sys.argv[1]
    print("IQ file=", path)
else:
    print("Specify IQ file.")
    exit()

with open(path, mode='rb') as f:
    x = np.fromfile(path, dtype=np.complex64)

### Sample rate and Symbol rate are required ###
### Please change them to match your system  ###
sample_rate = 1.024e6   # 1.024Mhz sampling
symbol_rate = 64e3      # 64kbps symbol rate
center_freq = 433.92e6  # 433.92MHz

# PSD
plt.psd(x, NFFT=1024, Fs=sample_rate/1e6, Fc=center_freq/1e6)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Relative power (dB)')
plt.title("Power Spectral Density")
plt.show()
"""
# Freq shift
N = len(x)
f_o = -57e3     # amount we need to shift by
t = np.arange(N)/sample_rate    # time vector
x = x * np.exp(2j*np.pi*f_o*t)  # down shift
"""

# Low-Pass Filter
### change it to match your system ###
freq_cutoff = 64e3
taps = firwin(numtaps=101, cutoff=freq_cutoff, fs=sample_rate)
x = np.convolve(x, taps, 'valid')

# Decimate by 2, now that we filtered and there wont be aliasing
x = x[::2]
sample_rate = 512e3 # 512ksample/s is integer ratio to symbol rate, 512/64=8

"""
### Sample rate and symbol rate are intger ratio
# Resample to 19kHz
x = resample_poly(x, 19, 25) # up, down
sample_rate = 19e3
"""

# Symbol sync, using what we did in sync chapter
samples = x # for the sake of matching the sync chapter
#
int_pol = 32
samples_interpolated = resample_poly(samples, int_pol, 1) # we'll use 32 as the interpolation factor, arbitrarily chosen, seems to work better than 16
# 512k/64k
sps = 8
mu = 0.01   # initial estimate of phase of sample
out = np.zeros(len(samples) + 10, dtype=np.complex64)
out_rail = np.zeros(len(samples) + 10, dtype=np.complex64) # stores values, each iteration we need the previous 2 values plus current value
i_in = 0    # input samples index
i_out = 2   # output index (let first two outputs be 0)
while i_out < len(samples) and i_in+int_pol < len(samples):
    out[i_out] = samples_interpolated[i_in*int_pol + int(mu*int_pol)] # grab what we think is the "best" sample
    out_rail[i_out] = int(np.real(out[i_out]) > 0) + 1j*int(np.imag(out[i_out]) > 0)
    x = (out_rail[i_out] - out_rail[i_out-2]) * np.conj(out[i_out-1])
    y = (out[i_out] - out[i_out-2]) * np.conj(out_rail[i_out-1])
    mm_val = np.real(y - x)
    mu += sps + 0.01*mm_val
    i_in += int(np.floor(mu))   # round down to nearest int since we are using it as an index
    mu = mu - np.floor(mu)      # remove the integer part of mu
    i_out += 1                  # increment output index
x = out[2:i_out]            # remove the first two, and anything after i_out (that was never filled out)

### After synbol sync, Circular Constellation is rotating because frequency offset remains,    
plt.scatter(np.real(x), np.imag(x), marker=".", alpha=0.25)
plt.title("After Symbol Sync")
plt.show()

# Fine freq sync =  Costas Loop
def cs_loop(x, alpha, beta):
    samples = x     # for the sake of matching the sync chapter
    N = len(samples)
    #print("alpha=", alpha, "beta=", beta)
    phase = 0
    freq = 0
    # These next two params is what to adjust, to make the feedback loop faster or slower (which impacts stability)
    #alpha = 8.0
    #beta = 0.002
    out = np.zeros(N, dtype=np.complex64)
    freq_log = []
    for i in range(N):
        out[i] = samples[i] * np.exp(-1j*phase) # adjust the input sample by the inverse of the estimated phase offset
        error = np.real(out[i]) * np.imag(out[i]) # This is the error formula for 2nd order Costas Loop (e.g. for BPSK)

        # Advance the loop (recalc phase and freq offset)
        freq += (beta * error)
        freq_log.append(freq * sample_rate / (2*np.pi)) # convert from angular velocity to Hz for logging
        phase += freq + (alpha * error)

        # Optional: Adjust phase so its always between 0 and 2pi, recall that phase wraps around every 2pi
        while phase >= 2*np.pi:
            phase -= 2*np.pi
        while phase < 0:
            phase += 2*np.pi
    x = out
    return(x, freq_log)

# Candidate values of alpha and beta
# Try another values depending on the result
alphas = [ 1, 2, 4, 8, 16 ]
betas = [ 0.002, 0.004, 0.008, 0.016, 0.032]

print("It takes a while, please be patient.")
### plot frequency log 
fig = plt.figure(tight_layout=True)
ax = fig.subplots(5,5)
for i,alpha in enumerate(alphas):
    for j,beta in enumerate(betas):
        xs = x.copy()
        _,freq_log = cs_loop(xs, alpha, beta)
        ax[i,j].plot(freq_log)
        ax[i,j].set_title("a="+str(alpha)+" b="+str(beta))
plt.show()
print("Remenber fast lock-in")

# plot constellation
fig = plt.figure(tight_layout=True)
ax = fig.subplots(5,5)
for i,alpha in enumerate(alphas):
    for j,beta in enumerate(betas):
        xs = x.copy()
        xr,fl = cs_loop(xs, alpha, beta)
        qty = len(xr) - 128     # temp lock-in symbols 
        ax[i,j].scatter(np.real(xr[-1*qty:]), np.imag(xr[ -1*qty:]), marker=".", alpha=0.25)
        ax[i,j].set_title("a="+str(alpha)+" b="+str(beta))
plt.show()
print("Remenber clear constellation")

### Set Best Value of alpha and beta
#alpha = 8
#beta = 0.016
alpha = float(input("Enter best alpha:"))
beta = float(input("Enter best beta:"))
xr,fl = cs_loop(x,alpha,beta)

#
plt.plot(fl[0:1000], label="freq_log, 0 to 1000 symbols")
plt.legend()
plt.show()
#
try:
    num_sym_lockin = int(input("Enter number of symbol lockin:"))
except:
    num_sym_lockin = 256

qty = len(xr) - num_sym_lockin 
plt.xlim([-1, 1])
plt.ylim([-0.2, 0.2])
plt.scatter(np.real(xr[-1*qty:]), np.imag(xr[ -1*qty:]), marker=".", alpha=0.25)
plt.title("a="+str(alpha)+" b="+str(beta)+" except first "+str(num_sym_lockin)+" symbols")
plt.show()

# Demod BPSK
bits = (np.real(x) > 0).astype(int) # 1's and 0's

# Differential decoding, so that it doesn't matter whether our BPSK was 180 degrees rotated without us realizing it
bits = (bits[1:] - bits[0:-1]) % 2
bits = bits.astype(np.uint8) # for decoder
print("bits:", bits)

# end of costas_loop_search.py