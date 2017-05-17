import fsk_lecim_phy as phyLayer
import fsk_lecim_mod as mod
import fsk_lecim_demod as demod

import matplotlib.pyplot as plt
import scipy.signal as sig
import numpy as np

if __name__ == '__main__':
	modulator = mod.fsk_lecim_modulator(sps=10, 
        modulationIndex=1.0, 
        Band169MHz=False, 
        phyLecimFskPreambleLength=4, 
        FCS=False, 
        dataWhitening=False, 
        pfsk=True, 
        phyPacketSize= 8)

	demodulator = demod.fsk_lecim_demodulator(sps=10, 
        modulationIndex=1.0, 
        Band169MHz=False, 
        FCS=True, 
        dataWhitening=False, 
        pfsk=True)

	mod_data = modulator.modulate_random(modulator.phyPacketSize)

	noise_power = 1
	bb = mod_data + (np.random.normal(size=len(mod_data), scale=np.sqrt(noise_power)) + 1j * np.random.normal(size=len(mod_data), scale=np.sqrt(noise_power)))

	demod_data = demodulator.demodulate(bb)

	bbrate = modulator.symbol_rate*modulator.sps

	#spectrum
	fig, ax = plt.subplots(3)
	f, Pwelch = sig.welch(mod_data, bbrate / 1000, nperseg=2048)
	ax[0].semilogy(np.fft.fftshift(f), np.fft.fftshift(Pwelch))
	ax[0].set_xlabel('Frequency / MHz')
	ax[0].set_ylabel('Spectrum')
	ax[0].set_title('Spectrum Baseband')
	#timesignal
	ax[1].plot(mod_data.real)
	ax[1].set_ylim([-1.2,1.2])
	ax[1].plot(mod_data.imag)
	ax[1].set_xlabel('time / '+ str(1.0/bbrate*1000.0) + 's')
	ax[1].set_ylabel('')
	ax[1].set_title('Timesignal')
	ax[2].plot(mod_data[0:10000].real, mod_data[0:10000].imag)
	ax[2].set_xlim([-2.2, 2.2])
	ax[2].set_ylim([-2.2, 2.2])
	ax[2].set_title('IQ')

	#spectrum with noise
	fig, ax = plt.subplots(1)
	f, Pwelch = sig.welch(bb, bbrate / 1000, nperseg=2048)
	ax.semilogy(np.fft.fftshift(f), np.fft.fftshift(Pwelch))
	ax.set_xlabel('Frequency / MHz')
	ax.set_ylabel('Spectrum')
	ax.set_title('Spectrum with noise')
	
	plt.show()


