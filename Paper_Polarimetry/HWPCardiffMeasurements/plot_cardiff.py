import numpy as np
import matplotlib.pyplot as plt

#dir = '/Users/macias/NIKA/Processing/Papers_Proposal_Notes/Paper_Polarimetry/HWPCardiffMeasurements/'
#dir = '/Users/ritacco/Nika/software/Processing/Papers_Proposal_Notes/Paper_Polarimetry/HWPCardiffMeasurements/'
dir = '/home/ritacco/nika/NIKA/Soft/Processing/Papers_Proposal_Notes/Paper_Polarimetry/HWPCardiffMeasurements/'
file = dir+'mesures_hwp_mesh_cardiff.txt'
data= np.loadtxt(file)



#NIKA BAND PASSES
import astropy.io.fits as fits
#dirbp = '/Users/macias/NIKA/Processing/Pipeline/Calibration/BP/'
dirbp = '/home/ritacco/nika/NIKA/Soft/Processing/Pipeline/Calibration/BP/'
filebp = 'NIKA_bandpass_Run8.fits'

hdu = fits.open(dirbp+filebp)
bp1 = hdu[1].data
bp2 = hdu[2].data

freq1 = bp1['FREQ']
trans1 = bp1['NIKATRANS'] 

freq2 = bp2['FREQ']
trans2 = bp2['NIKATRANS'] 

# FIGURES


fig = plt.figure(1)

plt.clf()

ax1 = fig.add_subplot(111)
ax1.plot(data[:,0],data[:,7],'ro-')
ax1.set_ylabel('$\Phi$[degrees]', size=26)
ax2 = ax1.twinx()
ax2.plot(freq1,trans1,'b')
ax2.plot(freq2,trans2,'c')
ax2.set_ylabel('Normalized Transmission',size=26)
ax1.set_xlabel('Frequency [GHz]',size=26)
ax1.set_ylim([30,205])
#plt.ylim([30,205])
plt.xlim([90,360])
ax1.tick_params(labelsize=26)
ax2.tick_params(labelsize=26)


plt.savefig(dir+'phase_shift_angle.pdf')
