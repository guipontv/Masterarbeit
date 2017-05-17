import numpy as np
import fsk_lecim_constants
import fsk_lecim_phy

from math import ceil, floor
from cmath import exp, pi

import commpy.channelcoding.convcode as cc

class fsk_lecim_demodulator(fsk_lecim_phy.physical_layer):
    def demodulate(self, data_in):
        if self.pfsk:
            PPDU = self.demodulator_pfsk(data_in)
        else:
            PPDU = self.demodulator_fsk(data_in)

        PPDU = self.PPDU_analyser(PPDU)
        PPDU[0] = self.deinterleaver(PPDU[0], True)
        self.PHR = self.fec_decoder(PPDU[0])
        print self.PHR
        self.PHR_analyser()
        
        PPDU[1] = self.deinterleaver(PPDU[1])
        PDU = self.fec_decoder(PPDU[1])
        PDU = self.zero_padding_remover(PDU)
        print PDU

    #demodulator FSK
    def demodulator_fsk(self, data_in):
        a = np.zeros((len(data_in),2), dtype = complex)
        data_out = np.zeros((int(len(data_in)/self.sps),), dtype = int)
        contribution = [0, 0, 0, 0]
        Z = [0, 0]

        for i in range(len(data_in)):
            a[i][0] = data_in[i]*exp(1j*2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))
            a[i][1] = data_in[i]*exp(1j*-2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))

        for k in range(int(len(data_out))):
            for p in range(self.sps):
                contribution[0] += (a[self.sps*k+p][0]).real 
                contribution[1] += (a[self.sps*k+p][0]).imag 
                contribution[2] += (a[self.sps*k+p][1]).real 
                contribution[3] += (a[self.sps*k+p][1]).imag 
            Z[0]= contribution[0]**2 + contribution[1]**2 #Z0
            Z[1]= contribution[2]**2 + contribution[3]**2 #Z1
            if Z[0] - Z[1] >= 0: #Z0-Z1 Threshold 0
                data_out[k] = 0
            else:
                data_out[k] = 1
            contribution = [0, 0, 0, 0]
            
        return data_out

    #Demodulator P-FSK
    def demodulator_pfsk(self, data_in):
        a = np.zeros((len(data_in),2), dtype = complex)
        data_out = np.zeros((int(len(data_in)/self.sps),), dtype = int)
        contribution = [0, 0, 0, 0]
        Z = [0, 0, 0, 0]
        delta = [0, 0]
        for i in range(len(data_in)):
            a[i][0] = data_in[i]*exp(1j*2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))
            a[i][1] = data_in[i]*exp(1j*-2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))
        for k in range(int(len(data_out))):
            for p in range(self.sps):
                contribution[0] += (a[self.sps*k+p][0]).real 
                contribution[1] += (a[self.sps*k+p][0]).imag 
                contribution[2] += (a[self.sps*k+p][1]).real 
                contribution[3] += (a[self.sps*k+p][1]).imag
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
    #deinterleave K    
    def deinterleave_k(self, k, Ndepth, lambda_):
        a = int((Ndepth-1-k)*lambda_-(Ndepth-1)*floor((Ndepth-1-k)*lambda_/float(Ndepth)))
        return a    

    #deinterleaver
    def deinterleaver(self, data_in, phr = False):
        data_out=np.zeros((len(data_in),),dtype = int)
        if(phr):
            for i in range(self.nphr):
                data_out[self.deinterleave_k(i,self.nphr,self.lambdaPhr)] = data_in[i]
        else:
            for m in range(self.nBlock):
                for k in range(self.npsdu):
                    data_out[m*self.npsdu+self.deinterleave_k(k,self.npsdu,self.lambdaPsdu)] = data_in[m*self.npsdu+k]
        return data_out

    #Viterbi/ FEC decoder
    def fec_decoder(self, data_in):
        memory = np.array([6])
        g_matrix = np.array([[91,121]]) # G(D) = [1+D^2+D^3+D^5+D^6, 1+D+D^2+D^3+D^6]
        trellis = cc.Trellis(memory, g_matrix)
        data_out = cc.viterbi_decode(data_in, trellis, tb_depth = int(len(data_in)/2))
        return data_out[:-6]

    #PPDU analyser 
    def PPDU_analyser(self, data_in):
        SHRlength = (int(floor(self.phyLecimFskPreambleLength))+3)*8
        return [data_in[SHRlength:SHRlength+44], data_in[SHRlength+44:]]

    #PHR analyser
    def PHR_analyser(self):
        pdu_len = 0
        if self.PHR[12]:
            self.FCS = True
        else:
            self.FCS = False
        if self.PHR[11]:
            self.dataWhitening = True
        else:
            self.dataWhitening = False
        for i in range(11):
            pdu_len = pdu_len + (self.PHR[i]<<i)
        parity_bit = self.PHR[0]
        for i in range(1,13):
            parity_bit = parity_bit^self.PHR[i]  
        if self.PHR[13]!=parity_bit:
            print 'Parity bit error'
        self.phyPacketSize = pdu_len
        self.nBlock = int(ceil((8*self.phyPacketSize+6)/(self.npsdu/2.0)))

    #zero padding remover
    def zero_padding_remover(self, data_in):
        npad = int((self.nBlock*(fsk_lecim_constants.Npsdu/2.0))-(8*self.phyPacketSize+6))
        self.npad = npad
        return data_in[:-npad]
