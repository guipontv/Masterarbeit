import fsk_lecim_phy as phyLayer
import fsk_lecim_mod as mod
import fsk_lecim_demod as demod
import numpy as np

if __name__ == '__main__':
	modulator = mod.fsk_lecim_modulator(sps=10, 
        modulationIndex=1.0, 
        Band169MHz=False, 
        phyLecimFskPreambleLength=4, 
        FCS=False, 
        dataWhitening=False, 
        pfsk=False, 
        phyPacketSize= 8)

	demodulator = demod.fsk_lecim_demodulator(sps=10, 
        modulationIndex=1.0, 
        Band169MHz=False, 
        FCS=False, 
        dataWhitening=False, 
        pfsk=False)

	mod_data = modulator.modulate_random(modulator.phyPacketSize)

	noise_power = 2
	bb = mod_data + (np.random.normal(size=len(mod_data), scale=np.sqrt(noise_power)) + 1j * np.random.normal(size=len(mod_data), scale=np.sqrt(noise_power)))

	demod_data = demodulator.demodulate(bb)