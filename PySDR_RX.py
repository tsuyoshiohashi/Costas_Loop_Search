#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
# pySDRbpsk64.py
#
# Reference : https://pysdr.org/content/rds.html
#
import numpy as np
from scipy.signal import resample_poly, firwin, bilinear, lfilter
import matplotlib.pyplot as plt

# Read in signal
x = np.fromfile('./43392bpsk64_1024MSample.iq', dtype=np.complex64)
sample_rate = 1.024e6
center_freq = 433.92e6

### FFT ###
N = len(x)
print(N)
dt = 1/sample_rate
F = np.fft.fft(x)
F = F/(N/2)
Amp = np.abs(F)
freq = np.fft.fftfreq(N,d=dt)
#print(freq)
plt.plot(freq, Amp, label="Amp1024k")
plt.legend()
plt.show()
###
"""
# Freq shift
N = len(x)
f_o = -57e3 # amount we need to shift by
t = np.arange(N)/sample_rate # time vector
x = x * np.exp(2j*np.pi*f_o*t) # down shift
"""

# Low-Pass Filter
taps = firwin(numtaps=101, cutoff=50e3, fs=sample_rate)
x = np.convolve(x, taps, 'valid')

# Decimate by 2, now that we filtered and there wont be aliasing
x = x[::2]
sample_rate = 512e3
"""
# Resample to 19kHz
x = resample_poly(x, 19, 25) # up, down
sample_rate = 19e3
"""
disp_fft = False
if(disp_fft):
    ### FFT ###
    N = len(x)
    print(N)
    dt = 1/sample_rate
    F = np.fft.fft(x)
    F = F/(N/2)
    Amp = np.abs(F)
    freq = np.fft.fftfreq(N,d=dt)
    #print(freq)
    plt.plot(freq, Amp, label="Amp512k")
    plt.legend()
    plt.show()
    ###

"""
###
fs = sample_rate
psd = np.fft.fftshift(np.abs(np.fft.fft(x)))
f = np.linspace(-fs/2.0, fs/2.0, len(psd))
plt.plot(f, psd, label="psd")
plt.legend()
plt.show()
###

max_freq = f[np.argmax(psd)]
print("max_freq=",max_freq)
Ts = dt # calc sample period
t = np.arange(0, Ts*N, Ts) # create time vector
x = x * np.exp(-1j*2*np.pi*max_freq*t/2.0)

F = np.fft.fft(x)
F = F/(N/2)
Amp = np.abs(F)
freq = np.fft.fftfreq(N,d=dt)
#print(freq)
plt.plot(freq, Amp, label="SFT")
plt.legend()
plt.show()
"""

# Symbol sync, using what we did in sync chapter
samples = x # for the sake of matching the sync chapter
#
int_pol = 32
samples_interpolated = resample_poly(samples, int_pol, 1) # we'll use 32 as the interpolation factor, arbitrarily chosen, seems to work better than 16
# 512k/64k
sps = 8
mu = 0.01 # initial estimate of phase of sample
out = np.zeros(len(samples) + 10, dtype=np.complex64)
out_rail = np.zeros(len(samples) + 10, dtype=np.complex64) # stores values, each iteration we need the previous 2 values plus current value
i_in = 0 # input samples index
i_out = 2 # output index (let first two outputs be 0)
while i_out < len(samples) and i_in+int_pol < len(samples):
    out[i_out] = samples_interpolated[i_in*int_pol + int(mu*int_pol)] # grab what we think is the "best" sample
    out_rail[i_out] = int(np.real(out[i_out]) > 0) + 1j*int(np.imag(out[i_out]) > 0)
    x = (out_rail[i_out] - out_rail[i_out-2]) * np.conj(out[i_out-1])
    y = (out[i_out] - out[i_out-2]) * np.conj(out_rail[i_out-1])
    mm_val = np.real(y - x)
    mu += sps + 0.01*mm_val
    i_in += int(np.floor(mu)) # round down to nearest int since we are using it as an index
    mu = mu - np.floor(mu) # remove the integer part of mu
    i_out += 1 # increment output index
x = out[2:i_out] # remove the first two, and anything after i_out (that was never filled out)

###
plt.scatter(np.real(x), np.imag(x))
plt.show()
###

def csloop(x, alpha, beta):
    # Fine freq sync
    samples = x # for the sake of matching the sync chapter
    N = len(samples)
    #print("alpha=", alpha, "beta=", beta)
    phase = 0
    freq = 0
    # These next two params is what to adjust, to make the feedback loop faster or slower (which impacts stability)
    #alpha = 32   #8.0
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
    #plt.plot(freq_log, '.-', label="freq_log")
    #plt.title("alpha="+str())
    #plt.legend()
    #plt.show()
#
alphas = [8, 16, 32, 64 ]
betas = [0.062, 0.125, 0.25, 0.5]
### plot frequency error 
fig = plt.figure(tight_layout=True)
ax = fig.subplots(4,4)
for i,alpha in enumerate(alphas):
    for j,beta in enumerate(betas):
        xs = x.copy()
        _,fl = csloop(xs, alpha, beta)
        ax[i,j].plot(fl)
        ax[i,j].set_title("a="+str(alpha)+" b="+str(beta))
        
        #qty = int(len(x) / 2) 
        #ax[i,j].scatter(np.real(x[-1*qty:]), np.imag(x[ -1*qty:]))
        #ax[i,j].set_title("a="+str(alpha)+" b="+str(beta))
plt.show()

### plot constalation
fig = plt.figure(tight_layout=True)
ax = fig.subplots(4,4)
for i,alpha in enumerate(alphas):
    for j,beta in enumerate(betas):
        xs = x.copy()
        xr,fl = csloop(xs, alpha, beta)
        #ax[i,j].plot(fl)
        #ax[i,j].set_title("a="+str(alpha)+" b="+str(beta))
        
        qty = int(len(xr) * 0.5) 
        ax[i,j].scatter(np.real(xr[-1*qty:]), np.imag(xr[ -1*qty:]))
        ax[i,j].set_title("a="+str(alpha)+" b="+str(beta))
plt.show()

####################
alpha = 8
beta = 0.062
xr,fl = csloop(x,alpha,beta)
#####################
plt.plot(fl, label="fl")
plt.legend()
plt.show()

qty = int(len(xr) * 0.5) 
plt.scatter(np.real(xr[-1*qty:]), np.imag(xr[ -1*qty:]))
plt.title("a="+str(alpha)+" b="+str(beta))
plt.show()
#####################
exit()
# Demod BPSK
bits = (np.real(x) > 0).astype(int) # 1's and 0's

# Differential decoding, so that it doesn't matter whether our BPSK was 180 degrees rotated without us realizing it
bits = (bits[1:] - bits[0:-1]) % 2
bits = bits.astype(np.uint8) # for decoder

# Constants
syndrome = [383, 14, 303, 663, 748]
offset_pos = [0, 1, 2, 3, 2]
offset_word = [252, 408, 360, 436, 848]

# see Annex B, page 64 of the standard
def calc_syndrome(x, mlen):
    reg = 0
    plen = 10
    for ii in range(mlen, 0, -1):
        reg = (reg << 1) | ((x >> (ii-1)) & 0x01)
        if (reg & (1 << plen)):
            reg = reg ^ 0x5B9
    for ii in range(plen, 0, -1):
        reg = reg << 1
        if (reg & (1 << plen)):
            reg = reg ^ 0x5B9
    return reg & ((1 << plen) - 1) # select the bottom plen bits of reg

# Initialize all the working vars we'll need during the loop
synced = False
presync = False

wrong_blocks_counter = 0
blocks_counter = 0
group_good_blocks_counter = 0

reg = np.uint32(0) # was unsigned long in C++ (64 bits) but numpy doesn't support bitwise ops of uint64, I don't think it gets that high anyway
lastseen_offset_counter = 0
lastseen_offset = 0

# the synchronization process is described in Annex C, page 66 of the standard */
bytes_out = []
for i in range(len(bits)):
    # in C++ reg doesn't get init so it will be random at first, for ours its 0s
    # It was also an unsigned long but never seemed to get anywhere near the max value
    # bits are either 0 or 1
    reg = np.bitwise_or(np.left_shift(reg, 1), bits[i]) # reg contains the last 26 rds bits. these are both bitwise ops
    if not synced:
        reg_syndrome = calc_syndrome(reg, 26)
        for j in range(5):
            if reg_syndrome == syndrome[j]:
                if not presync:
                    lastseen_offset = j
                    lastseen_offset_counter = i
                    presync = True
                else:
                    if offset_pos[lastseen_offset] >= offset_pos[j]:
                        block_distance = offset_pos[j] + 4 - offset_pos[lastseen_offset]
                    else:
                        block_distance = offset_pos[j] - offset_pos[lastseen_offset]
                    if (block_distance*26) != (i - lastseen_offset_counter):
                        presync = False
                    else:
                        print('Sync State Detected')
                        wrong_blocks_counter = 0
                        blocks_counter = 0
                        block_bit_counter = 0
                        block_number = (j + 1) % 4
                        group_assembly_started = False
                        synced = True
            break # syndrome found, no more cycles

    else: # SYNCED
        # wait until 26 bits enter the buffer */
        if block_bit_counter < 25:
            block_bit_counter += 1
        else:
            good_block = False
            dataword = (reg >> 10) & 0xffff
            block_calculated_crc = calc_syndrome(dataword, 16)
            checkword = reg & 0x3ff
            if block_number == 2: # manage special case of C or C' offset word
                block_received_crc = checkword ^ offset_word[block_number]
                if (block_received_crc == block_calculated_crc):
                    good_block = True
                else:
                    block_received_crc = checkword ^ offset_word[4]
                    if (block_received_crc == block_calculated_crc):
                        good_block = True
                    else:
                        wrong_blocks_counter += 1
                        good_block = False
            else:
                block_received_crc = checkword ^ offset_word[block_number] # bitwise xor
                if block_received_crc == block_calculated_crc:
                    good_block = True
                else:
                    wrong_blocks_counter += 1
                    good_block = False

            # Done checking CRC
            if block_number == 0 and good_block:
                group_assembly_started = True
                group_good_blocks_counter = 1
                bytes = bytearray(8) # 8 bytes filled with 0s
            if group_assembly_started:
                if not good_block:
                    group_assembly_started = False
                else:
                    # raw data bytes, as received from RDS. 8 info bytes, followed by 4 RDS offset chars: ABCD/ABcD/EEEE (in US) which we leave out here
                    # RDS information words
                    # block_number is either 0,1,2,3 so this is how we fill out the 8 bytes
                    bytes[block_number*2] = (dataword >> 8) & 255
                    bytes[block_number*2+1] = dataword & 255
                    group_good_blocks_counter += 1
                    #print('group_good_blocks_counter:', group_good_blocks_counter)
                if group_good_blocks_counter == 5:
                    #print(bytes)
                    bytes_out.append(bytes) # list of len-8 lists of bytes
            block_bit_counter = 0
            block_number = (block_number + 1) % 4
            blocks_counter += 1
            if blocks_counter == 50:
                if wrong_blocks_counter > 35: # This many wrong blocks must mean we lost sync
                    print("Lost Sync (Got ", wrong_blocks_counter, " bad blocks on ", blocks_counter, " total)")
                    synced = False
                    presync = False
                else:
                    print("Still Sync-ed (Got ", wrong_blocks_counter, " bad blocks on ", blocks_counter, " total)")
                blocks_counter = 0
                wrong_blocks_counter = 0

# Constants
syndrome = [383, 14, 303, 663, 748]
offset_pos = [0, 1, 2, 3, 2]
offset_word = [252, 408, 360, 436, 848]

# see Annex B, page 64 of the standard
def calc_syndrome(x, mlen):
    reg = 0
    plen = 10
    for ii in range(mlen, 0, -1):
        reg = (reg << 1) | ((x >> (ii-1)) & 0x01)
        if (reg & (1 << plen)):
            reg = reg ^ 0x5B9
    for ii in range(plen, 0, -1):
        reg = reg << 1
        if (reg & (1 << plen)):
            reg = reg ^ 0x5B9
    return reg & ((1 << plen) - 1) # select the bottom plen bits of reg

# Initialize all the working vars we'll need during the loop
synced = False
presync = False

wrong_blocks_counter = 0
blocks_counter = 0
group_good_blocks_counter = 0

reg = np.uint32(0) # was unsigned long in C++ (64 bits) but numpy doesn't support bitwise ops of uint64, I don't think it gets that high anyway
lastseen_offset_counter = 0
lastseen_offset = 0

# the synchronization process is described in Annex C, page 66 of the standard */
bytes_out = []
for i in range(len(bits)):
    # in C++ reg doesn't get init so it will be random at first, for ours its 0s
    # It was also an unsigned long but never seemed to get anywhere near the max value
    # bits are either 0 or 1
    reg = np.bitwise_or(np.left_shift(reg, 1), bits[i]) # reg contains the last 26 rds bits. these are both bitwise ops
    if not synced:
        reg_syndrome = calc_syndrome(reg, 26)
        for j in range(5):
            if reg_syndrome == syndrome[j]:
                if not presync:
                    lastseen_offset = j
                    lastseen_offset_counter = i
                    presync = True
                else:
                    if offset_pos[lastseen_offset] >= offset_pos[j]:
                        block_distance = offset_pos[j] + 4 - offset_pos[lastseen_offset]
                    else:
                        block_distance = offset_pos[j] - offset_pos[lastseen_offset]
                    if (block_distance*26) != (i - lastseen_offset_counter):
                        presync = False
                    else:
                        print('Sync State Detected')
                        wrong_blocks_counter = 0
                        blocks_counter = 0
                        block_bit_counter = 0
                        block_number = (j + 1) % 4
                        group_assembly_started = False
                        synced = True
            break # syndrome found, no more cycles

    else: # SYNCED
        # wait until 26 bits enter the buffer */
        if block_bit_counter < 25:
            block_bit_counter += 1
        else:
            good_block = False
            dataword = (reg >> 10) & 0xffff
            block_calculated_crc = calc_syndrome(dataword, 16)
            checkword = reg & 0x3ff
            if block_number == 2: # manage special case of C or C' offset word
                block_received_crc = checkword ^ offset_word[block_number]
                if (block_received_crc == block_calculated_crc):
                    good_block = True
                else:
                    block_received_crc = checkword ^ offset_word[4]
                    if (block_received_crc == block_calculated_crc):
                        good_block = True
                    else:
                        wrong_blocks_counter += 1
                        good_block = False
            else:
                block_received_crc = checkword ^ offset_word[block_number] # bitwise xor
                if block_received_crc == block_calculated_crc:
                    good_block = True
                else:
                    wrong_blocks_counter += 1
                    good_block = False

            # Done checking CRC
            if block_number == 0 and good_block:
                group_assembly_started = True
                group_good_blocks_counter = 1
                bytes = bytearray(8) # 8 bytes filled with 0s
            if group_assembly_started:
                if not good_block:
                    group_assembly_started = False
                else:
                    # raw data bytes, as received from RDS. 8 info bytes, followed by 4 RDS offset chars: ABCD/ABcD/EEEE (in US) which we leave out here
                    # RDS information words
                    # block_number is either 0,1,2,3 so this is how we fill out the 8 bytes
                    bytes[block_number*2] = (dataword >> 8) & 255
                    bytes[block_number*2+1] = dataword & 255
                    group_good_blocks_counter += 1
                    #print('group_good_blocks_counter:', group_good_blocks_counter)
                if group_good_blocks_counter == 5:
                    #print(bytes)
                    bytes_out.append(bytes) # list of len-8 lists of bytes
            block_bit_counter = 0
            block_number = (block_number + 1) % 4
            blocks_counter += 1
            if blocks_counter == 50:
                if wrong_blocks_counter > 35: # This many wrong blocks must mean we lost sync
                    print("Lost Sync (Got ", wrong_blocks_counter, " bad blocks on ", blocks_counter, " total)")
                    synced = False
                    presync = False
                else:
                    print("Still Sync-ed (Got ", wrong_blocks_counter, " bad blocks on ", blocks_counter, " total)")
                blocks_counter = 0
                wrong_blocks_counter = 0
                