import fsk_lecim_phy as phyLayer
import fsk_lecim_mod as mod
import fsk_lecim_demod as demod

import matplotlib.pyplot as plt
import scipy.signal as sig
import numpy as np
from cmath import exp, pi

#Compare data
def assertData(data1, data2):
	if len(data1) != len(data2):
		print "Length error " + str(len(data1)) + " != " + str(len(data2))
		return -2
	ctr = 0
	for i in range(len(data1)):
		if data1[i]!=data2[i]:
			ctr += 1
			print "Bit error " + str(i) + " position"
	if ctr != 0:
		print str(ctr) + "bit(s) error"
		return -1
	print 'Equals'
	return 0

if __name__ == '__main__':
	pfsk = True
	index = 1.0
	modulator = mod.fsk_lecim_modulator(sps=20, 
        modulationIndex=index, 
        Band169MHz=False, 
        phyLecimFskPreambleLength=4, 
        FCS=False, 
        dataWhitening=False, 
        pfsk=pfsk, 
        phyPacketSize=32)

	demodulator = demod.fsk_lecim_demodulator(sps=20, 
        modulationIndex=index, 
        Band169MHz=False, 
        FCS=False, 
        dataWhitening=False, 
        pfsk=pfsk)

	mod_data = modulator.modulate_random(modulator.phyPacketSize)

	fading = exp(-1j*3*pi/2)
	noise_power = 1
	bb = mod_data[1] * fading + (np.random.normal(size=len(mod_data[1]), scale=np.sqrt(noise_power)) + 1j * np.random.normal(size=len(mod_data[1]), scale=np.sqrt(noise_power)))
	demod_data = demodulator.demodulate(bb,modulator.phyPacketSize)
	demod_data1 = demodulator.demodulate_coherent(bb,modulator.phyPacketSize)

	#assert DATA
	print "asserting PHR"
	print modulator.PHR
	print 'non coherent'
	print demod_data[0]
	test1 = assertData(modulator.PHR, demod_data[0])
	print 'coherent'
	print demod_data1[0]
	test2 = assertData(modulator.PHR, demod_data1[0])
	
	print "asserting PDU"
	print mod_data[0]
	print 'non coherent'
	print demod_data[1]
	test3 = assertData(mod_data[0], demod_data[1])
	print 'coherent'
	print demod_data1[1]
	test4 = assertData(mod_data[0], demod_data1[1])

	bbrate = modulator.symbol_rate*modulator.sps

	#spectrum no noise
	fig, ax = plt.subplots(3)
	f, Pwelch = sig.welch(mod_data[1], bbrate / 1000, nperseg=2048)
	ax[0].semilogy(np.fft.fftshift(f), np.fft.fftshift(Pwelch))
	ax[0].set_xlabel('Frequency / MHz')
	ax[0].set_ylabel('Spectrum')
	ax[0].set_title('Spectrum Baseband')

	#time signal
	ax[1].plot(mod_data[1][0:240].real)
	ax[1].set_ylim([-1.2,1.2])
	ax[1].plot(mod_data[1][0:240].imag)
	ax[1].set_xlabel('time / '+ str(1.0/bbrate*1000.0) + 's')
	ax[1].set_ylabel('')
	ax[1].set_title('\nTimesignal')
	ax[2].plot(mod_data[1][0:10000].real, mod_data[1][0:10000].imag)
	ax[2].set_xlim([-2.2, 2.2])
	ax[2].set_ylim([-2.2, 2.2])
	ax[2].set_title('\nIQ')

	#spectrum with noise
	fig, ax = plt.subplots(1)
	f, Pwelch = sig.welch(bb, bbrate / 1000, nperseg=2048)
	ax.semilogy(np.fft.fftshift(f), np.fft.fftshift(Pwelch))
	ax.set_xlabel('Frequency / MHz')
	ax.set_ylabel('Spectrum')
	ax.set_title('Spectrum with noise')


	plt.show()
