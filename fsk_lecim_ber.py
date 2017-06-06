import fsk_lecim_phy as phyLayer
import fsk_lecim_mod as mod
import fsk_lecim_demod as demod

import matplotlib.pyplot as plt
import scipy.signal as sig
import scipy.special as sp
import numpy as np
from cmath import exp, pi


def BER(modulator, demodulator):
	mod_data = modulator.modulate_random(modulator.phyPacketSize)
	SNR = []
	ctrnoncoherent = []
	ctrcoherent = []
	theoryBer = []
	ctr0 = 0
	ctr1 = 0
	fading_factor = 1
	for i in range(8):
		noise_power = 10**(-i/10.0)
		SNR.append(i)
		for j in range(2):
			print 'step ' + str(i) + ' ' + str(j)
			channel = modulator.rayleigh_fading(mod_data[1], fading_factor, noise_power)
			bb = mod_data[1] * channel[0] + channel[1]			
			demod_data = demodulator.demodulate(bb, modulator.phyPacketSize)
			channel = modulator.rayleigh_fading(mod_data[1], fading_factor, noise_power)
			bb = mod_data[1] * channel[0] + channel[1]			
			demod_data1 = demodulator.demodulate_coherent(bb, modulator.phyPacketSize)
			ctr0 += sum(abs(mod_data[0]-demod_data[1]))
			ctr1 += sum(abs(mod_data[0]-demod_data1[1]))
		ctrcoherent.append(ctr0/float(len(2*mod_data[0])))
		ctrnoncoherent.append(ctr1/float(2*len(mod_data[0])))
		theoryBer.append(0.5*sp.erfc(np.sqrt(10**(i/10.0))))
		print noise_power
		print SNR[-1]
		print theoryBer[-1]
		print ctrcoherent[-1]
		print ctrnoncoherent[-1]
		ctr0 = 0
		ctr1 = 0
		if ctrcoherent[-1] == 0 and ctrnoncoherent[-1] == 0:
			break
	fig, ax = plt.subplots(1)
	ax.semilogy(SNR, ctrcoherent)
	ax.semilogy(SNR, ctrnoncoherent)
	ax.semilogy(SNR, theoryBer)
	ax.set_xlabel('Eb/N0 (dB)')
	ax.set_ylabel('Bit error')
	ax.set_title('BER')
	plt.show()

if __name__ == '__main__':
	pfsk = False
	index = 1.0
	modulator = mod.fsk_lecim_modulator(sps=20, 
        modulationIndex=index, 
        Band169MHz=False, 
        phyLecimFskPreambleLength=4, 
        FCS=False, 
        dataWhitening=False, 
        pfsk=pfsk, 
        phyPacketSize=2047)

	demodulator = demod.fsk_lecim_demodulator(sps=20, 
        modulationIndex=index, 
        Band169MHz=False, 
        FCS=False, 
        dataWhitening=False, 
        pfsk=pfsk,
        phyPacketSize=2047)

	BER(modulator, demodulator)