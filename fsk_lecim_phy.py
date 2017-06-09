import fsk_lecim_constants
import numpy as np

from math import ceil

class physical_layer:
	def __init__(self, 
		sps=10, 
		modulationIndex=1.0, 
		Band169MHz=False, 
		phyLecimFskPreambleLength=4, 
		FCS=False, 
		dataWhitening=False, 
		pfsk=False, 
		phyLecimFskSpreading = False, phyLecimFskSpreadingAlternating = False, phyLecimFskSpreadingFactor = 2, 
		phyPacketSize=1):
	
		#SHR parameter
		self.preamble = fsk_lecim_constants.preamble	
		self.SFD = fsk_lecim_constants.SFD
		#PHR parameter	
		self.phyLecimFskPreambleLength = phyLecimFskPreambleLength
		self.FCS = FCS
		self.dataWhitening = dataWhitening
		self.phyPacketSize = phyPacketSize if phyPacketSize <= fsk_lecim_constants.aMaxPhyPacketSize else fsk_lecim_constants.aMaxPhyPacketSize
		#SHR & PHR generation
		self.SHR = self.gen_SHR()
		self.PHR = self.gen_PHR()
		#interleaver parameter
		self.nPhr = fsk_lecim_constants.nPhr
		self.nPsdu = fsk_lecim_constants.nPsdu
		self.lambdaPhr = fsk_lecim_constants.lambdaPhr
		self.lambdaPsdu = fsk_lecim_constants.lambdaPsdu
		self.nBlock = int(ceil((8*self.phyPacketSize+6)/(self.nPsdu/2.0)))
		#spreading
		self.phyLecimFskSpreading = phyLecimFskSpreading
		self.phyLecimFskSpreadingFactor = phyLecimFskSpreadingFactor
		self.phyLecimFskSpreadingAlternating = phyLecimFskSpreadingAlternating
		#zero padding
		self.nPad = int((self.nBlock*(fsk_lecim_constants.nPsdu/2))-(8*phyPacketSize+6))
		#modulation
		self.pfsk = pfsk
		self.sps  = sps
		self.modulationIndex = modulationIndex
		self.Band169MHz = Band169MHz
		self.symbol_rate = self.symbol_rate()
		self.freq_dev = self.freq_dev()
	
	#Compute symbol rate
	def symbol_rate(self):
		if self.Band169MHz:
			if self.modulationIndex == 1:
				R = 12500 #bit/s
				print "Symbol rate is 12.5 kb/s"
			elif self.modulationIndex == 0.5:
				R = 25000
				print "Symbol rate is 25 kb/s"
			else:
				raise Exception("Modulation index must be 0.5 or 1.0 for the 169 MHz band")
		else:
			if self.modulationIndex == 0.5:
				R = 37500
				print "Symbol rate is 37.5 kb/s"
			elif self.modulationIndex == 1.0:
				R = 25000
				print "Symbol rate is 25 kb/s"
			elif self.modulationIndex == 2.0:
				R = 12500
				print "Symbol rate is 12.5 kb/s"
			else:
				raise Exception("Modulation index must be 0.5 or 1.0 or 2.0 for bands > 169 MHz band")
		return R		

	#Frequency deviation
	def freq_dev(self):
		fdev = (self.modulationIndex*self.symbol_rate)/2
		print "Frequency deviation is " + str(fdev) + " Hz"
		return fdev

	#SHR generator
	def gen_SHR(self):
		if self.phyLecimFskPreambleLength < 4 and self.phyLecimFskPreambleLength > 64:
			raise Exception("Preamble length must be between 4 and 64")
		SHR = np.tile(self.preamble,self.phyLecimFskPreambleLength)
		SHR = np.concatenate((SHR,self.SFD))
		return SHR

	#PHR generator
	def gen_PHR(self):
		PHR = np.zeros((16, ), dtype=int)
		if self.FCS:
			PHR[12] = 1
		else:
			PHR[12] = 0
		if self.dataWhitening:
			PHR[11] = 1
		else:
			PHR[11] = 0
		payl_len_bitstring = '{0:011b}'.format(self.phyPacketSize)
		payl_len_list = [int(payl_len_bitstring[i],2) for i in range(0,len(payl_len_bitstring))]
		PHR[0:11] = np.flipud(payl_len_list)
		parity_bit = PHR[0]
		for i in range(1,13):
			parity_bit = parity_bit^PHR[i]
		PHR[13] = parity_bit
		return PHR