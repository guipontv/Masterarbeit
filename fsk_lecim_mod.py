import numpy as np
import fsk_lecim_constants
import fsk_lecim_phy

from math import floor
from cmath import exp, pi

class fsk_lecim_modulator(fsk_lecim_phy.physical_layer):
    def modulate_random(self, n):
        self.phyPacketSize = n
        self.PHR = self.gen_PHR()
        PDU = np.random.randint(0, 2, size=(n * 8,))
        return self.modulate(PDU)

    def modulate(self, data_in):
        print data_in
        print self.PHR
        PHR = self.fec_encoder(self.PHR)
        PHR = self.interleaver(PHR, True)
        
        PDU = self.zero_padding(data_in)
        PDU = self.fec_encoder(PDU)
        PDU = self.interleaver( PDU)
        PSDU = self.mux(self.SHR, PHR, PDU)
        PSDU = self.mapper(PSDU)
        PSDU = self.modulator(PSDU)
        return PSDU

    #PDU length in bytes
    def pdu_len(self, data_in):
        pdu_len = len(data_in)
        if((pdu_len%8)!= 0):
            print 'PDU length is not a multiple of 8, padding zero to compensate'
            data_in = np.concatenate((data_in, np.zeros((8-(pdu_len%8), ),dtype = int)))
        return len(data_in)/8

    #zero padding
    def zero_padding(self, data_in):
        npad = self.nPad
        data_out = []
        data_out = np.concatenate((data_in, np.zeros((npad, ),dtype = int)))
        return data_out

    #FEC
    def fec_encoder(self, data_in):
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
    def interleave_k(self, k, Ndepth, lambda_):
        a = int((((Ndepth-1-k)%lambda_)*Ndepth/float(lambda_))+floor((Ndepth-1-k)/(float(lambda_))))
        return a    

    #interleaver
    def interleaver(self, data_in, phr = False):
        data_out=np.zeros((len(data_in),),dtype = int)
        if(phr):
            for i in range(self.nphr):
                data_out[self.interleave_k(i,self.nphr,self.lambdaPhr)] = data_in[i]
        else:
            for m in range(self.nBlock):
                for k in range(self.npsdu):
                    data_out[m*self.npsdu+self.interleave_k(k,self.npsdu,self.lambdaPsdu)] = data_in[m*self.npsdu+k]
        return data_out

    #MUX
    def mux(self, shr, phr, psdu):

        return np.concatenate((shr, phr, psdu))
 
    #mapper
    def mapper(self, data_in):
        data_out = np.zeros((len(data_in),),dtype = int)
        if self.pfsk:
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
    def modulator(self, data_in):
        data_out = np.zeros((len(data_in)*self.sps,), dtype = complex)
        for i in range(len(data_in)):
            for k in range(self.sps):
                data_out[self.sps*i+k] = abs(data_in[i])*exp(1j*2*pi*data_in[i]*self.freq_dev*(self.sps*i+k)/(self.sps*self.symbol_rate))
        return data_out