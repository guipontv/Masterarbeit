import numpy as np

aMaxPhyPacketSize = 2047 # The maximum PSDU size (in octets) the PHY shall be able to receive, gnuradio doesn't accept more than 16383 imput data
preamble = np.array([0,1,0,1,0,1,0,1], dtype=int)
SFD = np.array([0,1,1,1,0,0,0,0,1,1,1,0,1,1,1,0,1,1,0,1,0,0,1,0], dtype=int)
lambdaPhr = 4
lambdaPsdu = 6
nPhr = 44 # nDepth (nPhr = 4*11)
nPsdu = 72 # nDepth (nPsdu = 6*12)
spreadingAlternating_0= np.array([0,1], dtype=int)
spreadingAlternating_1= np.array([1,0], dtype=int)
spreadingNonAlternating8_0 = np.array([1,0,1,1,0,0,0,1], dtype=int)
spreadingNonAlternating8_1 = np.array([0,1,0,0,1,1,1,0], dtype=int)
spreadingNonAlternating16_0 = np.array([0,0,1,0, 0,0,1,1, 1,1,0,1, 0,1,1,0], dtype=int)
spreadingNonAlternating16_1 = np.array([1,1,0,1, 1,1,0,0, 0,0,1,0, 1,0,0,1], dtype=int)