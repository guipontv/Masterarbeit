import fsk_lecim_phy
import numpy as np
from math import ceil, floor
from cmath import exp, pi

#PDU length in bytes
def pdu_len_block(data_in):
	pdu_len = len(data_in)
	if((pdu_len%8)!= 0):
		print 'PDU length is not a multiple of 8, padding zero to compensate'
		data_in = np.concatenate((data_in, np.zeros((8-(pdu_len%8), ),dtype = int)))
	return len(data_in)/8

#zero padding
def zero_padding_block(phy, data_in):
	npad = phy.nPad
	data_out = []
	data_out = np.concatenate((data_in, np.zeros((npad, ),dtype = int)))
	return data_out

#FEC
def fec_block(data_in):
	poly = np.zeros((6,),dtype=int)
	data_in = np.concatenate((data_in, np.zeros((6,),dtype = int)))
	data_out = np.zeros((2*len(data_in)),dtype = int)
	for i in range(len(data_in)):
		data_out[2*i] = data_in[i]^poly[1]^poly[2]^poly[4]^poly[5]
		data_out[2*i+1] = data_in[i]^poly[0]^poly[1]^poly[2]^poly[5]
		poly = np.insert(poly,0, data_in[i])
		poly = np.delete(poly,6)
	return data_out

#interleave K    
def interleave_k(k, Ndepth, lambda_):
    a = int((((Ndepth-1-k)%lambda_)*Ndepth/float(lambda_))+floor((Ndepth-1-k)/(float(lambda_))))
    return a    

#interleaver
def interleaver_block(phy,data_in,phr = False):
	data_out=np.zeros((len(data_in),),dtype = int)
	if(phr):
		for i in range(phy.nphr):
			data_out[interleave_k(i,phy.nphr,phy.lambdaPhr)] = data_in[i]
	else:
		for m in range(phy.nBlock):
			for k in range(phy.npsdu):
				data_out[m*phy.npsdu+interleave_k(k,phy.npsdu,phy.lambdaPsdu)] = data_in[m*phy.npsdu+k]
	return data_out

#deinterleave K    
def deinterleave_k(k, Ndepth, lambda_):
    a = int((Ndepth-1-k)*lambda_-(Ndepth-1)*floor((Ndepth-1-k)*lambda_/float(Ndepth)))
    return a    

#deinterleaver
def deinterleaver_block(phy, data_in,phr = False):
	data_out=np.zeros((len(data_in),),dtype = int)
	if(phr):
		for i in range(phy.nphr):
			data_out[deinterleave_k(i,phy.nphr,phy.lambdaPhr)] = data_in[i]
	else:
		for m in range(phy.nBlock):
			for k in range(phy.npsdu):
				data_out[m*phy.npsdu+deinterleave_k(k,phy.npsdu,phy.lambdaPsdu)] = data_in[m*phy.npsdu+k]
	return data_out

#MUX
def mux_block(shr, phr, psdu):
	return np.concatenate((shr, phr, psdu))

#mapper
def mapper_block(phy, data_in):
	data_out = np.zeros((len(data_in),),dtype = int)
	if phy.pfsk:
		for i in range(int(floor(len(data_in)/2.0))):
			if data_in[2*i] == 0 and data_in[2*i+1] == 0:
				data_out[2*i] = -1
				data_out[2*i+1]  = 0
			if data_in[2*i] == 0 and data_in[2*i+1] == 1:
				data_out[2*i] = 0
				data_out[2*i+1] = -1
			if data_in[2*i] == 1 and data_in[2*i+1] == 0:
				data_out[2*i] = 1
				data_out[2*i+1] = 0
			if data_in[2*i] == 1 and data_in[2*i+1] == 1:
				data_out[2*i] = 0
				data_out[2*i+1] = 1
	else:
		for i in range(len(data_in)):
			if data_in[i]:
				data_out[i] = 1
			else:
				data_out[i] = -1
	return data_out

#modulator
def modulator_block(phy, data_in):
	data_out = np.zeros((len(data_in)*phy.sps,), dtype = complex)
	for i in range(len(data_in)):
		for k in range(phy.sps):
			data_out[phy.sps*i+k] = abs(data_in[i])*exp(1j*2*pi*data_in[i]*phy.freq_dev*(phy.sps*i+k)/(phy.sps*phy.symbol_rate))
	return data_out

#demodulator FSK
def demodulator_fsk_block(phy, data_in):
	a = np.zeros((len(data_in),2), dtype = complex)
	data_out = np.zeros((int(len(data_in)/phy.sps),), dtype = int)
	contribution = [0, 0, 0, 0]
	Z = [0, 0]

	for i in range(len(data_in)):
		a[i][0] = data_in[i]*exp(1j*2*pi*phy.freq_dev*i/(phy.sps*phy.symbol_rate))
		a[i][1] = data_in[i]*exp(1j*-2*pi*phy.freq_dev*i/(phy.sps*phy.symbol_rate))

	for k in range(int(len(data_out))):
		for p in range(phy.sps):
			contribution[0] += (a[phy.sps*k+p][0]).real 
			contribution[1] += (a[phy.sps*k+p][0]).imag 
			contribution[2] += (a[phy.sps*k+p][1]).real 
			contribution[3] += (a[phy.sps*k+p][1]).imag 
		Z[0]= contribution[0]**2 + contribution[1]**2 #Z0
		Z[1]= contribution[2]**2 + contribution[3]**2 #Z1
		if Z[0] - Z[1] >= 0: #Z0-Z1 Threshold 0
			data_out[k] = 0
		else:
			data_out[k] = 1
		contribution = [0, 0, 0, 0]
		
	return data_out

#Demodulator P-FSK
def demodulator_pfsk_block(phy, data_in):
	a = np.zeros((len(data_in),2), dtype = complex)
	data_out = np.zeros((int(len(data_in)/phy.sps),), dtype = int)
	contribution = [0, 0, 0, 0]
	Z = [0, 0, 0, 0]
	delta = [0, 0]
	for i in range(len(data_in)):
		a[i][0] = data_in[i]*exp(1j*2*pi*phy.freq_dev*i/(phy.sps*phy.symbol_rate))
		a[i][1] = data_in[i]*exp(1j*-2*pi*phy.freq_dev*i/(phy.sps*phy.symbol_rate))
	for k in range(int(len(data_out))):
		for p in range(phy.sps):
			contribution[0] += (a[phy.sps*k+p][0]).real 
			contribution[1] += (a[phy.sps*k+p][0]).imag 
			contribution[2] += (a[phy.sps*k+p][1]).real 
			contribution[3] += (a[phy.sps*k+p][1]).imag
		if (k%2) == 0:
			Z[0]= contribution[0]**2 + contribution[1]**2 #Z0(2k)
			Z[1]= contribution[2]**2 + contribution[3]**2 #Z1(2k)
		else:
			Z[2]= contribution[0]**2 + contribution[1]**2 #Z0(2k+1)
			Z[3]= contribution[2]**2 + contribution[3]**2 #Z1(2k+1)
			
			delta = [Z[0]+Z[1], Z[2]+Z[3]]

			if delta[0] - delta[1] >= 0: #Position bit
				data_out[k] = 0
				if Z[0]-Z[1]>=0:
					data_out[k-1] = 0
				else:
					data_out[k-1] = 1
			else:
				data_out[k] = 1
				if Z[2]-Z[3]>=0:
					data_out[k-1] = 0
				else:
					data_out[k-1] = 1
		contribution = [0, 0, 0, 0]
	return data_out

#PPDU analyser : preamble length, SHR, PHR
def PPDU_analyser(phy, data_in):
	preambleLength = 0
	for i in range(int(len(data_in)/2.0)):
		if data_in[2*i] == 0 and data_in[2*i+1] == 1:
			preambleLength += 0.25
		else:
			break
	SHRlength = (int(floor(preambleLength))+3)*8
	phy.phyLecimFskPreambleLength = int(floor(preambleLength))
	phy.SHR = data_in[:SHRlength]
	phy.PHR = data_in[SHRlength:SHRlength+44]
	return data_in[SHRlength+44:]

#PHR analyser
def PHR_analyser(phy,phr):
	pdu_len = 0
	if phr[12]:
		phy.FCS = True
	else:
		phy.FCS = False
	if phr[11]:
		phy.dataWhitening = True
	else:
		phy.dataWhitening = False
	for i in range(11):
		pdu_len = pdu_len + (phr[i]<<i)
	phy.phyPacketSize = pdu_len
	phy.PHR = phy.gen_PHR()
	if phy.PHR[13]!=phr[13]:
		print 'Parity bit error'
	phy.nBlock = int(ceil((8*phy.phyPacketSize+6)/(phy.npsdu/2.0)))


	



