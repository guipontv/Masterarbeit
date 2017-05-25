import fsk_lecim_phy as phyLayer
import fsk_lecim_mod as mod
import fsk_lecim_demod as demod

import matplotlib.pyplot as plt
import scipy.signal as sig
import numpy as np
from cmath import exp, pi


def BER(modulator, demodulator):
	mod_data = modulator.modulate_random(modulator.phyPacketSize)
	SNR = []
	ctrnoncoherent = []
	ctrcoherent = []
	ctr = 0
	for i in range(3):
		print 'step ' + str(i)
		fading = 1
		noise_power = 5.0 - 0.2*i 
		print noise_power
		SNR.append(10*np.log(5/noise_power))
		bb = mod_data[1] * fading + (np.random.normal(size=len(mod_data[1]), scale=np.sqrt(noise_power)) + 1j * np.random.normal(size=len(mod_data[1]), scale=np.sqrt(noise_power)))
		demod_data = demodulator.demodulate(bb, modulator.phyPacketSize)
		bb = mod_data[1] * fading + (np.random.normal(size=len(mod_data[1]), scale=np.sqrt(noise_power)) + 1j * np.random.normal(size=len(mod_data[1]), scale=np.sqrt(noise_power)))
		demod_data1 = demodulator.demodulate_coherent(bb, modulator.phyPacketSize)
		ctr = sum(abs(mod_data[0]-demod_data[1]))/float(len(mod_data[0]))
		print ctr 
		ctrcoherent.append(ctr)
		ctr = sum(abs(mod_data[0]-demod_data1[1]))/float(len(mod_data[0]))
		print ctr
		ctrnoncoherent.append(ctr)

	print ctrcoherent
	print ctrnoncoherent
	print SNR
	"""	fig, ax = plt.subplots(2)
	ax[0].semilogy(ctrcoherent)
	ax[0].set_xlabel('SNR')
	ax[0].set_ylabel('BER coherent')
	ax[0].set_title('BER')
	ax[1].semilogy(ctrnoncoherent)
	ax[1].set_xlabel('SNR')
	ax[1].set_ylabel('BER non coherent ')
	ax[1].set_title('BER')
	plt.show()"""

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
        phyPacketSize=1024)

	demodulator = demod.fsk_lecim_demodulator(sps=20, 
        modulationIndex=index, 
        Band169MHz=False, 
        FCS=False, 
        dataWhitening=False, 
        pfsk=pfsk,
        phyPacketSize=1024)

	BER(modulator, demodulator)

	


