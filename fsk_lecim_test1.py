import fsk_lecim_phy as phyLayer
import fsk_lecim_blocks as bl
import numpy as np

#put result in a string
def write_result(data_in):
	data_txt = ""
	for elt in data_in:
		a = str(elt)
		a = a[1:-1]
		data_txt = data_txt + a + "\n"
	return data_txt

#put data in a string
def write_data(data_in):
    data_txt = "["
    for elt in data_in:
        a = str(elt)
        data_txt = data_txt + a + ", "
    data_txt = data_txt[:-2] + "]\n\n"
    return data_txt

#Compare data
def assertData(data1, data2):
	if len(data1) != len(data2):
		print "Size error"
		return 1
	for i in range(len(data1)):
		if data1[i]!=data2[i]:
			print str(i) + "th bit error"
			return 1
	print 'Equality'
	return 0

if __name__ == '__main__':
    """
    TRANSMITTER TEST
    """
    print 'Test transmitter chain :'
    print 'Test pdu_len_block'
    PDU = np.random.randint(2, size=256)
    PDUinit = PDU
    pdu_len = bl.pdu_len_block(PDU)
    print pdu_len

    print 'Test physical-layer'
    phy = phyLayer.physical_layer( 
        sps=10, 
        modulationIndex=1.0, 
        Band169MHz=False, 
        phyLecimFskPreambleLength=4, 
        FCS=False, 
        dataWhitening=False, 
        pfsk=False, 
        phyPacketSize= pdu_len)
    SHR = phy.SHR
    PHR = phy.PHR
    print 'SHR:'
    print phy.SHR
    print 'PHR:'
    print phy.PHR

#PHR processing
    print 'data processing PHR'

    PHR = bl.fec_block(PHR)

    print 'PHR after FEC'
    print PHR
    PHRafterfec = PHR
    PHR = bl.interleaver_block(phy, PHR, True)

    print 'PHR after interleaver'
    print PHR
    PHRafterinterleaver = PHR

#PDU processing    
    print 'data processing PDU'
    PDU = bl.zero_padding_block(phy,PDU)

    print 'PDU after zero padding'
    PDUafterzeropadding = PDU

    PDU = bl.fec_block(PDU)

    print 'PDU after FEC'
    PDUafterfec = PDU

    PDU = bl.interleaver_block(phy, PDU)

    print 'PDU after interleaver'
    PDUafterinterleaver = PDU

#MUX
    PPDUbeforemapping = bl.mux_block(SHR,PHR,PDU)

#Mapping
    PPDU = bl.mapper_block(phy, PPDUbeforemapping)
    print 'PPDU after mapping'
    PPDUaftermapping = PPDU
    
#Modulation
    PPDU = bl.modulator_block(phy,PPDU)

    print 'PPDU modulated'
    print 'len modulated PPDU'
    print len(PPDU)

    noise_power = 1   
    modulatedPPDU = PPDU + (np.random.normal(size=len(PPDU), scale=np.sqrt(noise_power)) + 1j * np.random.normal(size=len(PPDU), scale=np.sqrt(noise_power)) )

    """
    RECEIVER TEST
    """
    print 'Test receiver chain :'
    phyReceiver = phyLayer.physical_layer(sps=10, modulationIndex=1.0, Band169MHz=False, pfsk=False)

#Demodulation
    if phy.pfsk:
    	print 'PFSK'
    	PPDU = bl.demodulator_pfsk_block(phy,PPDU)
    else:
    	print 'FSK'
    	PPDU = bl.demodulator_fsk_block(phy,PPDU)
    demodulatedPPDU = PPDU

    print 'Demodulated PPDU test'
    assertData(PPDU,PPDUbeforemapping)
    print len(PPDU)

    print 'PPDU analyser'
    PDU = bl.PPDU_analyser(phyReceiver, PPDU)
    print len(PDU)
    extractedPDU = PDU

#PHR Processing
#Deinterleaving
    phyReceiver.PHR = bl.deinterleaver_block(phyReceiver,phyReceiver.PHR, True)
    print 'PHR deinterleaving'
    extractedPHR = phyReceiver.PHR
    assertData(phyReceiver.PHR, PHRafterfec)
#Viterbi 
    phyReceiver.PHR = bl.fec_decoder_block(phyReceiver.PHR)
    decodedPHR = phyReceiver.PHR
    print "decoded PHR"
    print decodedPHR

    print "PHR"
    print phy.PHR

    assertData(phyReceiver.PHR, phy.PHR)

#PHR extraction

    bl.PHR_analyser(phyReceiver,phy.PHR)

#PDU processing
#Deinteleaving
    PDU = bl.deinterleaver_block(phyReceiver,PDU, False)
    print 'PDU deinterleaving'
    assertData(PDU, PDUafterfec)
#Viterbi

    decodedPDU = bl.fec_decoder_block(PDU)
    
    print "decoded PDU"
    print decodedPDU

    print "PDUafterzeropadding"
    print PDUafterzeropadding

    assertData(PDUafterzeropadding, decodedPDU)

    PDU = bl.zero_padding_remover_block(phyReceiver,decodedPDU)
    assertData(PDU, PDUinit)
#Write file result    
    result = open("result.txt","w")
    data_txt = write_result(modulatedPPDU)
    result.write(data_txt)
    result.close()
#Write file data 
    datafile = open("data.txt","w")

    datafile.write("SHR\n")
    data_txt = write_data(phy.SHR)
    datafile.write(data_txt)

    datafile.write("PHR\n")
    data_txt = write_data(phy.PHR)
    datafile.write(data_txt)

    datafile.write("PHR after FEC\n")
    data_txt = write_data(PHRafterfec)
    datafile.write(data_txt)

    datafile.write("PHR after interleaver\n")
    data_txt = write_data(PHRafterinterleaver)
    datafile.write(data_txt)

    datafile.write("PDU\n")
    data_txt = write_data(PDUinit)
    datafile.write(data_txt)

    datafile.write("PDU after zero padding\n")
    data_txt = write_data(PDUafterzeropadding)
    datafile.write(data_txt)

    datafile.write("PDU after FEC\n")
    data_txt = write_data(PDUafterfec)
    datafile.write(data_txt)

    datafile.write("PDU after interleaver\n")
    data_txt = write_data(PDUafterinterleaver)
    datafile.write(data_txt)

    datafile.write("PPDU\n")
    data_txt = write_data(PPDUbeforemapping)
    datafile.write(data_txt)

    datafile.write("PPDU after mapping\n")
    data_txt = write_data(PPDUaftermapping)
    datafile.write(data_txt)
    """
    datafile.write("Modulated PPDU\n")
    data_txt = write_data(modulatedPPDU)
    datafile.write(data_txt)
    """
    datafile.write("Demodulated PPDU\n")
    data_txt = write_data(demodulatedPPDU)
    datafile.write(data_txt)

    datafile.write("Extracted SHR\n")
    data_txt = write_data(phyReceiver.SHR)
    datafile.write(data_txt)

    datafile.write("Extracted PHR\n")
    data_txt = write_data(extractedPHR)
    datafile.write(data_txt)

    datafile.write("Extracted PDU\n")
    data_txt = write_data(extractedPDU)
    datafile.write(data_txt)

    datafile.write("Deinterleaved PHR\n")
    data_txt = write_data(phyReceiver.PHR)
    datafile.write(data_txt)

    datafile.write("Deinterleaved PDU\n")
    data_txt = write_data(PDU)
    datafile.write(data_txt)

    datafile.close()



