import numpy as np
import fsk_lecim_constants
import fsk_lecim_phy

from math import ceil, floor, cos
from cmath import exp, pi, phase
from scipy import signal
import commpy.channelcoding.convcode as cc

class fsk_lecim_demodulator(fsk_lecim_phy.physical_layer):
    def demodulate(self, data_in, phyPacketSize):
        if self.pfsk:
            PPDU = self.demodulator_pfsk(data_in)
        else:
            PPDU = self.demodulator_fsk(data_in)

        PPDU = self.PPDU_analyser(PPDU)
        PPDU[0] = self.deinterleaver(PPDU[0], True)
        PHR = self.fec_decoder(PPDU[0])
        self.PHR = PHR
        self.PHR_analyser()
        if self.phyPacketSize != phyPacketSize:
            print 'good PDU length'
            self.phyPacketSize = phyPacketSize
            self.nBlock = int(ceil((8*self.phyPacketSize+6)/(self.npsdu/2.0)))
        PPDU[1] = self.deinterleaver(PPDU[1])
        PDU = self.fec_decoder(PPDU[1])
        PDU = self.zero_padding_remover(PDU)
        return [PHR, PDU]

    def demodulate_coherent(self, data_in, phyPacketSize):
        if self.pfsk:
            PPDU = self.demodulator_pfsk_coherent(data_in)
        else:
            PPDU = self.demodulator_fsk_coherent(data_in)

        PPDU = self.PPDU_analyser(PPDU)
        PPDU[0] = self.deinterleaver(PPDU[0], True)
        PHR = self.fec_decoder(PPDU[0])
        self.PHR = PHR
        self.PHR_analyser()
        if self.phyPacketSize != phyPacketSize:
            print 'good PDU length'
            self.phyPacketSize = phyPacketSize
            self.nBlock = int(ceil((8*self.phyPacketSize+6)/(self.npsdu/2.0)))
        PPDU[1] = self.deinterleaver(PPDU[1])
        PDU = self.fec_decoder(PPDU[1])
        PDU = self.zero_padding_remover(PDU)
        return [PHR, PDU]

    #demodulator FSK correlator (non coherent)
    def demodulator_fsk(self, data_in):
        a = np.zeros((len(data_in),2), dtype = complex)
        data_out = np.zeros((int(len(data_in)/self.sps),), dtype = int)
        sum0 = [0, 0]
        Z = [0, 0]
        for i in range(len(data_in)):
            a[i][0] = data_in[i]*exp(1j*2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))
            a[i][1] = data_in[i]*exp(1j*-2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))
        for k in range(int(len(data_out))):
            for p in range(self.sps):
                sum0[0] += a[self.sps*k+p][0]
                sum0[1] += a[self.sps*k+p][1]
            Z[0]= abs(sum0[0])**2 #Z0
            Z[1]= abs(sum0[1])**2 #Z1
            if Z[0] - Z[1] >= 0: #Z0-Z1 Threshold 0
                data_out[k] = 0
            else:
                data_out[k] = 1
            sum0 = [0, 0]
        return data_out

    #demodulator FSK correlator (coherent)
    def demodulator_fsk_coherent(self, data_in):
        a = np.zeros((len(data_in),2), dtype = complex)
        data_out = np.zeros((int(len(data_in)/self.sps),), dtype = int)
        sum0 = [0, 0]
        Z = [0, 0]
        d = [0, 0]
        for i in range(len(data_in)):
            a[i][0] = data_in[i]*exp(1j*2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))
            a[i][1] = data_in[i]*exp(1j*-2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))
            
        for k in range(int(len(data_out))):
            for p in range(self.sps):
                sum0[0] += a[self.sps*k+p][0]
                sum0[1] += a[self.sps*k+p][1]
            phioff = phase(sum0[0]*d[0]+sum0[1]*d[1])
            sum0[0] = sum0[0] * exp(-1j*phioff)
            sum0[1] = sum0[1] * exp(-1j*phioff)
            Z[0]= sum0[0].real #Z0
            Z[1]= sum0[1].real #Z1
            if Z[0] - Z[1] >= 0: #Z0-Z1 Threshold 0
                data_out[k] = 0
                d[0] = self.sps
            else:
                data_out[k] = 1
                d[1] = self.sps
            sum0 = [0, 0]
        return data_out

    #Demodulator P-FSK correlator (non coherent)
    def demodulator_pfsk(self, data_in):
        a = np.zeros((len(data_in),2), dtype = complex)
        data_out = np.zeros((int(len(data_in)/self.sps),), dtype = int)
        sum0 = [0, 0]
        Z = [0, 0, 0, 0]
        delta = [0, 0]
        for i in range(len(data_in)):
            a[i][0] = data_in[i]*exp(1j*2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))
            a[i][1] = data_in[i]*exp(1j*-2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))
        for k in range(int(len(data_out))):
            for p in range(self.sps):
                sum0[0] += a[self.sps*k+p][0] 
                sum0[1] += a[self.sps*k+p][1]
            if (k%2) == 0:
                Z[0]= abs(sum0[0])**2 #Z0(2k)
                Z[1]= abs(sum0[1])**2 #Z1(2k)
            else:
                Z[2]= abs(sum0[0])**2 #Z0(2k+1)
                Z[3]= abs(sum0[1])**2 #Z1(2k+1)
                
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
            sum0 = [0, 0]
        return data_out

    #Demodulator P-FSK correlator (coherent)
    def demodulator_pfsk_coherent(self, data_in):
        a = np.zeros((len(data_in),2), dtype = complex)
        data_out = np.zeros((int(len(data_in)/self.sps),), dtype = int)
        sum0 = [0, 0]
        Z = [0, 0, 0, 0]
        d = [0, 0]
        delta = [0, 0]
        for i in range(len(data_in)):
            a[i][0] = data_in[i]*exp(1j*2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))
            a[i][1] = data_in[i]*exp(1j*-2*pi*self.freq_dev*i/(self.sps*self.symbol_rate))
        for k in range(int(len(data_out))):
            for p in range(self.sps):
                sum0[0] += a[self.sps*k+p][0] 
                sum0[1] += a[self.sps*k+p][1]
            phioff = phase(sum0[0]*d[0]+sum0[1]*d[1])
            sum0[0] = sum0[0] * exp(-1j*phioff)
            sum0[1] = sum0[1] * exp(-1j*phioff)
            
            if (k%2) == 0:
                Z[0]= sum0[0].real #Z0(2k)
                Z[1]= sum0[1].real #Z1(2k)
            else:
                Z[2]= sum0[0].real #Z0(2k+1)
                Z[3]= sum0[1].real #Z1(2k+1)
                
                delta = [Z[0]+Z[1], Z[2]+Z[3]]

                if delta[0] - delta[1] >= 0: #Position bit
                    data_out[k] = 0
                    if Z[0]-Z[1]>=0:
                        data_out[k-1] = 0
                        d[0] = self.sps
                    else:
                        data_out[k-1] = 1
                        d[1] = self.sps
                else:
                    data_out[k] = 1
                    if Z[2]-Z[3]>=0:
                        data_out[k-1] = 0
                        d[0] = self.sps
                    else:
                        data_out[k-1] = 1
                        d[1] = self.sps
            sum0 = [0, 0]
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
    #PLL 
    def phase_loop(self, data_in):
        phi_ = np.zeros((len(data_in),), dtype = complex)
        phi = 0
       
        wn = 0.11
        zeta = 0.707
        K = 1000

        t1 = K/(wn*wn)
        t2 = 2*zeta/wn
        b0 = (4*K/t1)*(1+0.5*t2)
        b1 = 8*K/t1
        b2 = (4*K/t1)*(1-0.5*t2)
        a1 = -2.0
        a2 =  1.0

        v0 = 0
        v1 = 0
        v2 = 0
        
        diff = []
        for i in range(len(data_in)):
            delta = phase(data_in[i]*np.conj(exp(1j*phi)))
            diff.append(delta)
            #print str(i) +' Phi estimate ' + str(phi) + ' Error ' + str(delta)
            v2 = v1
            v1 = v0
            v0 = delta - v1*a1 - v2*a2
            phi = v0*b0 + v1*b1 + v2*b2
            phi_[i] = exp(1j*phi)

        return [phi_, diff]
