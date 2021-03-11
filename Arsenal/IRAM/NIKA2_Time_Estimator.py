#!/usr/bin/env python
##############################################################################################
### Created by: P. Garcia                                                                  ###
### Revised by: A. Ritacco                                                                 ###
### Revised by: B. Ladjelate                                                               ###
### Title: NIKA2_Time_Estimator.py                                                         ###
### Purpose: Calculation of the total integration time for observing proposals             ###
### Creation Date: 2016.JULY.28                                                            ###
### Date of last revision: 2018.JAN.20                                                     ###
### Revision history:                                                                      ###
###   - 2016.JULY.29: comments from C. Kramer & A. Sievers implemented                     ###
###   - 2016.AUG.03 : changed formula to N. Billot expression                              ###
###   - 2017.JAN.29 : - fpix set to 0.75 (commisi. results)                                ###
###                   - NEFD0 @ 2 mm set to 15 mJy/sqrt(hz) (commisi. results)             ###
###                   - overhead set as used defined parameter between 1.5 - 2.0           ###
###   - 2017.JULY.28: - According to the NIKA2 technical paper submitted on July, 3rd 2017 ###
###                   - fpix @1mm, 2mm set to 0.84 and 0.90, respectively                  ###
###                   - NEFD0 @1mm, 2mm set to 20 and 6 mJy/sqrt(hz), respectively         ###
###                   - FWHM set to 11.2 @1mm and 17.7 @2mm                                ###
###   - 2017.SEPT.28  - An error in the data software led sensitivities too optimistics    ###
###                   - NEFD0 @1mm, 2mm set to 40 and 10 mJy/sqrt(hz), respectively        ###
###   - 2018.JAN.20:  - NEFD0 @1mm, 2mm set to 33 and 8 mJy/sqrt(hz), respectively         ###
###   - 2019.FEB.18:  - Corrected a mistake in the calculation of the mapping speed        ###
###   - 2020.JUL.22:  - Implemented a warning regarding the scan size                      ###
##############################################################################################
import numpy as np, sys
import math
import sys
import os
os.system("rm -f *.py~")

#version = 'v 2017.JAN.29'
#version = 'v 2017.JULY.28'
#version = 'v 2017.SEPT.28'
#version = 'v 2018.JAN.20'
#version = 'v 2019.FEB.18'
version = 'v 2020.JUL.22'


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#print ''
print bcolors.UNDERLINE + bcolors.HEADER + bcolors.BOLD + 'Time-estimator-NIKA2 '+version+' '+bcolors.ENDC
#print ''
print bcolors.HEADER + bcolors.BOLD + ' See the "Guidelines for observing time estimates with the NIKA2 continuum camera' + bcolors.ENDC
print bcolors.HEADER + bcolors.BOLD + ' at the IRAM-30m Telescope" for details on used parameters and calculations.\n' + bcolors.ENDC

################################
##### CONVERSION FACTORS #######
################################

##########  K_cmb ############
def b_v_cmb(freq,tcmb):
    h    = 6.626070040*10**(-34) # (J.s) Planck constant
    k    = 1.38064852*10**(-23)  # (J/K) Boltzmann constant
    c    = 299792458             # (m/s^2) light speed
    freq = freq*10**9            # Hz    
    bv = 1.0/(((2*h*freq**3)/((c**2)*(np.e**((h*freq)/(k*tcmb)) - 1)))*(np.e**((h*freq)/(k*tcmb))/(np.e**((h*freq)/(k*tcmb)) - 1))*(h*freq/(k*tcmb**2))) / 10**20
    return bv

##########  Ysz ##############
def y_sz(freq,tcmb):
    h    = 6.626070040*10**(-34) # (J.s) Planck constant
    k    = 1.38064852*10**(-23)  # (J/K) Boltzmann constant
    c    = 299792458             # (m/s^2) light speed
    freq = freq*10**9            # Hz    
    Ysz  = 1.0
    bv   = 1.0/(((2*h*freq**3)/((c**2)*(np.e**((h*freq)/(k*tcmb)) - 1)))*(np.e**((h*freq)/(k*tcmb))/(np.e**((h*freq)/(k*tcmb)) - 1))*(h*freq/(k*tcmb**2)))
    bv   = 1.0/bv
    y_sz = (bv*tcmb)*(((h*freq/(k*tcmb))*((np.e**((h*freq)/(k*tcmb)) + 1)/(np.e**((h*freq)/(k*tcmb)) -1))) -4)*Ysz /10**(-20)
    y_sz = 1.0/y_sz
    return y_sz

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%% FIXED GENERAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

################################################
##### BEAM SIZES AND CONVERSION TO RADIANS #####
################################################
res_band1_spec        = 11.2  #From NIKA2 1st paper #MEMO # 12.0  # [arcsec]
res_band2_spec        = 17.7  #From NIKA2 1st paper #MEMO # 20.0  # [arcsec]
teta_band1_spec       = res_band1_spec*math.pi/180.0/3600.0 # in sr units 
teta_band2_spec       = res_band2_spec*math.pi/180.0/3600.0 # in sr units 
area_teta_band1_spec  = ((teta_band1_spec/2.0)**2)*math.pi  # in sr units 
area_teta_band2_spec  = ((teta_band2_spec/2.0)**2)*math.pi  # in sr units 
############################
##### OTHER PARAMETERS #####
############################
#fpix_spec             = 0.75
fpix_spec_1mm         = 0.84    # As presented at the NIKA2 consortium meeting
fpix_spec_2mm         = 0.90    # As presented at the NIKA2 consortium meeting
bv_1                  = 0.075   # band1, do not change
cv_1                  = 0.001   # band1, do not change
bv_2                  = 0.025   # band2, do not change
cv_2                  = 0.001   # band2, do not change
tiempo                = 1.0     # [hours]
FoV                   = 6.5     # [arcmin] diameter
FoVArea               = ((FoV/2.0)**2)*math.pi # [arcmin^2]
narrays1              = 2
narrays2              = 1
verbose_flag          = 0
#foverhead             = 2.0 # value fixed for all proposals in 2016
##########################################
##### NEFDo VALUES FROM OBSERVATIONS #####
##########################################
NEFD_band1_spec       = 33 #Changed in Jan. 20 #From NIKA2 1st paper arXiv:1707.00908  #MEMO 
NEFD_band2_spec       = 8  #Changed in Jan. 20 #From NIKA2 1st paper arXiv:1707.00908  #comissioning results 2017

#######################################
##### SAMPLING ANGULAR VELOCITIES #####
#######################################
sampl_vel_fast = 48.0    # [arcsec/second] 
sampl_vel_mid  = 24.0    # [arcsec/second] 
sampl_vel_slow = 12.0    # [arcsec/second] 
sampl_rate     = 23.84   # [hz]
sampl_rate_pol = 2*23.84 # [hz]
slewing_loss   = 1.0     # to account for slewing of single OTF line
##########################################
##########################################
##########################################

rango = len(sys.argv)
for i in range(rango):

    if sys.argv[i] == ("--help") or rango == 1:

        print ' USAGE:\n'
        print '   python NIKA2_Time_Estimator.py --help '
        print '   python NIKA2_Time_Estimator.py --band 1 --rms 2.00 --pwv 2 --elevation 40 --Xsize 6.5 --Ysize 6.5 --filter 1.0 --overhead 1.5'
        print '   python NIKA2_Time_Estimator.py --band 2 --rms 1.00 --pwv 4 --elevation 50 --Xsize 15  --Ysize 15  --filter 2.0 --overhead 2.0 --verbose '
        print '\n'
        print ' OPTIONS:\n'
        print '          help =>  This help.'
        print '          band =>  Set 1 or 2 for the 1 mm or the 2 mm bands, respectively.'
        print '           rms =>  Wanted flux density per beam. Any value above the confusion limit in [mJy/beam].'
        print '           pwv =>  Precipitable water vapor in [mm].'
        print '     elevation =>  Values from 15 to < 83 [deg].'
        print '   Xsize Ysize =>  Map lengths. Xsize and Ysize are in [arcmin]. Minimum map size is 6.5x6.5 [arcmin^2] for homogeneous RMS noise distribution.'
        print '        filter =>  Factor for post-processing noise filtering. Values from 1.0 (point-like source) to 2.0 (extended bright emission).'
        print '      overhead =>  Factor for telescope overheads between 1.5 and 2.0 (slewing, pointing, focusing, calibration), i.e. all telescope time which is not spend integrating on-source.'
        print '       verbose =>  Set to get list of parameters used in the calculations, RMS noise unit conversion, and allowed OTF scan speeds.\n'
        sys.exit()

    if sys.argv[i] == ("--band"):
        band  = int(sys.argv[i+1]) 
        if band == 1:
            narrays = narrays1
            mili    = 1.2
        if band == 2:
            narrays = narrays2
            mili    = 2.0
        if band < 1 or band > 2:
            print '    %%%%%%%%%%%%%%%%%%%%%%%'
            print bcolors.BOLD +  bcolors.FAIL +'    Band %2i is not defined.' % (band) + bcolors.ENDC
            print '    %%%%%%%%%%%%%%%%%%%%%%%'
            print ''
            sys.exit()

    if sys.argv[i] == ("--rms"):
        rms = float(sys.argv[i+1])

    if sys.argv[i] == ("--pwv"):
        pwv = float(sys.argv[i+1])

    if sys.argv[i] == ("--elevation"):
        elevation = float(sys.argv[i+1])
        if elevation < 15 or elevation > 83:
            print '    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            print bcolors.BOLD +  bcolors.FAIL +'    Elevation %3i [degrees] is outside the telescope`s limits: 15 - 83 [degrees].' % (elevation) + bcolors.ENDC
            print '    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            print ''
            sys.exit()

    if sys.argv[i] == ("--Xsize"):
        dx = float(sys.argv[i+1])

    if sys.argv[i] == ("--Ysize"):
        dy = float(sys.argv[i+1])

    if sys.argv[i] == ("--filter"):
        filtering = float(sys.argv[i+1])
        if filtering < 1 or filtering > 2:
            print '    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            print bcolors.BOLD +  bcolors.FAIL +'    Filter value %2.1f is outside the standard limits: 1.0 - 2.0' % (filtering) + bcolors.ENDC
            print '    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            print ''
            sys.exit()

    if sys.argv[i] == ("--overhead"):
        foverhead = float(sys.argv[i+1])
        if foverhead < 1.5 or foverhead > 2:
            print '    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            print bcolors.BOLD +  bcolors.FAIL +'    Overhead value %2.1f is outside the standard limits: 1.5 - 2.0' % (foverhead) + bcolors.ENDC
            print '    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            print ''
            sys.exit()

    if sys.argv[i] == ("--verbose"):
        verbose_flag = 1

area_obs        = dx*dy

if area_obs < 4.0:
    print '    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print bcolors.BOLD +  bcolors.FAIL +'    Introduced map size is %3.1f [arcmin^2]. Minimum value is 4.0 [arcmin^2]' % (area_obs) + bcolors.ENDC
    print '    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print ''
    sys.exit()

###############################################################################################
### SIGMA ESTIMATION CALCULATION, EFFECTIVE NEFD FROM PWV GIVEN VALUES AND UNITS CONVERSION ###
###############################################################################################
        
if band == 1: 
    tau = bv_1*pwv + cv_1
    NEFD_spec_eff      = NEFD_band1_spec*(np.e**(tau/np.sin(elevation*math.pi/180.)))
    rms_point_MJy_spec = rms/area_teta_band1_spec/10**9
    rms_point_cmb_spec = rms_point_MJy_spec*b_v_cmb(260,2.726)*10**6
    rms_point_ys_spec  = rms_point_MJy_spec*y_sz(260,2.726)*10**6
    rms_point_ys_spec  = np.abs(rms_point_ys_spec)
    #########################################
    ### ratio of areas with fraction of valid pixel in the 2mm band ###       
    #########################################
    factor_area_spec   = (1+(area_obs/(FoVArea*fpix_spec_1mm))) # from N. Billot Doc.
    ##########################################################
    ### Integration Time calculations for given conditions ###
    ##########################################################
    t_spec     = ((NEFD_spec_eff*filtering/rms)**2)*factor_area_spec*(foverhead)/3600.0 #hours
    s_map_spec = ((FoVArea*fpix_spec_1mm)/(((NEFD_spec_eff*filtering)**2)*foverhead))*3600.0            # arcmin^2 / hour^-1/ mJy^-2 
if band == 2:
    tau = bv_2*pwv + cv_2
    NEFD_spec_eff      = NEFD_band2_spec*(np.e**(tau/np.sin(elevation*math.pi/180.)))
    rms_point_MJy_spec = rms/area_teta_band2_spec/10**9
    rms_point_cmb_spec = rms_point_MJy_spec*b_v_cmb(150,2.726)*10**6   
    rms_point_ys_spec  = rms_point_MJy_spec*y_sz(150,2.726)*10**6     
    rms_point_ys_spec  = np.abs(rms_point_ys_spec)
    #########################################
    ### ratio of areas with fraction of valid pixel in the 2mm band ###       
    #########################################
    factor_area_spec   = (1+(area_obs/(FoVArea*fpix_spec_2mm))) # from N. Billot Doc.
    ##########################################################
    ### Integration Time calculations for given conditions ###
    ##########################################################
    t_spec     = ((NEFD_spec_eff*filtering/rms)**2)*factor_area_spec*(foverhead)/3600.0 #hours
    s_map_spec = ((FoVArea*fpix_spec_2mm)/(((NEFD_spec_eff*filtering)**2)*foverhead))*3600.0            # arcmin^2 / hour^-1/ mJy^-2 

#########################################
#### TIME PER OTF LINE TO REACH TOTAL ###
#########################################

T_SLOW = slewing_loss*((dx*60.0)/sampl_vel_slow)/60.0  # [minutes] 
T_MID  = slewing_loss*((dx*60.0)/sampl_vel_mid)/60.0   # [minutes] 
T_FAST = slewing_loss*((dx*60.0)/sampl_vel_fast)/60.0  # [minutes] 

if verbose_flag == 1:

    print ''
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print '%%%%%%%%%% General Parameters Used for the Calculations %%%%%%%%%%%'
    print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print ''
    print '------------------------------------------------------------------------'
    print '|      Parameters      |    Band 1 (1.2 mm)    |    Band 2 (2.0 mm)    |'
    print '------------------------------------------------------------------------'
    print '|      opacity         |  %6.3fxpwv + %6.3f  |  %6.3fxpwv + %6.3f  |' % (bv_1,cv_1,bv_2,cv_2)
    print '------------------------------------------------------------------------'
    print '| HPBW [arcsec]        |       %6.1f          |       %6.1f          |' % (res_band1_spec,res_band2_spec)
    print '------------------------------------------------------------------------'
    print '| NEFDo [mJy.s^0.5]    |       %6.1f          |       %6.1f          |' % (NEFD_band1_spec,NEFD_band2_spec) 
    print '------------------------------------------------------------------------'
    print '| fpix                 |       %6.2f          |       %6.2f          |' % (fpix_spec_1mm, fpix_spec_2mm)
    print '------------------------------------------------------------------------'
    if band  == 1:
        print '| rms [mJy/beam]       |       %6.2f          |                       |' % (rms)
    if band  == 2:
        print '| rms [mJy/beam]       |                       |       %6.2f          |' % (rms)
    print '------------------------------------------------------------------------'
    print '| FoV [arcmin]         |                   %6.1f                      |' % (FoV) 
    print '------------------------------------------------------------------------'
    print '| h-filtering          |                   %6.2f                      |' % (filtering) 
    print '------------------------------------------------------------------------'
    print '| h-overhead           |                   %6.2f                      |' % (foverhead) 
    print '------------------------------------------------------------------------'
    print '| OTF slow [arcsec/s]  |                    %4i                       |' % (sampl_vel_slow)
    print '------------------------------------------------------------------------'
    print '| OTF mid  [arcsec/s]  |                    %4i                       |' % (sampl_vel_mid)
    print '------------------------------------------------------------------------'
    print '| OTF fast [arcsec/s]  |                    %4i                       |' % (sampl_vel_fast)
    print '------------------------------------------------------------------------'
    print '| Dump [hz]            |                   %6.2f                      |' % (sampl_rate)
    print '------------------------------------------------------------------------'
    print '| Dump POL [hz]        |                   %6.2f                      |' % (sampl_rate_pol)
    print '------------------------------------------------------------------------'
    print ''

    print '     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print '     %%%%%%%%%%%%%%%%%   Units Conversion   %%%%%%%%%%%%%%%%%'
    print '     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print ''    
    print '     --------------------------------------------------------'
    print '     |         (User Defined)             |    Band %1i       |' % (band)
    print '     --------------------------------------------------------'
    print '     |rms [mJy/beam]                      |  %8.2f       |' % (rms) 
    print '     --------------------------------------------------------'
    print '     |rms point-like source [MJy/sr]      |     %6.2f      |'  % (rms_point_MJy_spec) 
    print '     |rms point-like source [K_CMB]x10^-6 |     %5i       |'  % (rms_point_cmb_spec) 
    print '     |*rms point.like source [Ysz]x10^-6  |     %5i       |'  % (rms_point_ys_spec) 
    print '     --------------------------------------------------------'
    print '     |* for Ysz = 1.0                                       |'
    print '     --------------------------------------------------------'
    print ''


#############################
#### TABLE FINAL RESULTS ####
#############################
        
print '     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'    
print '     %%%%%%%%%%%%%%%%%%%%%%   Results  %%%%%%%%%%%%%%%%%%%%%%'    
print '     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'    
print ''    
print '     for: El = %2i (deg), PWV = %1i (mm), tau = %4.2f, rms = %6.2f'  %  (elevation,pwv,tau,rms) 
print '     Area = %4.1f (arcmin^2), Filter = %4.1f, Overhead = %4.1f'  %  (area_obs,filtering,foverhead) 
print '  ---------------------------------------------------------------'
print '  |    (User Defined)      |                          |         |'
print '  ---------------------------------------------------------------'
print '  | Band                   |                          |    %1i    |'  % (band) 
print '  ---------------------------------------------------------------'
if band  == 1:
    print '  | Total Integration Time at 1mm |         [hours]          | %6.1f  |'  % (t_spec)
    print '  ---------------------------------------------------------------'
    print '  | Mapping speed          | [arcmin^-2/hour/mJy^-2 ] | %6.1f  |'  % (s_map_spec) 
if band  == 2:
    print '  | Total Integration Time at 2mm |         [hours]          | %6.1f  |'  % (t_spec)
    print '  ---------------------------------------------------------------'
    print '  | Mapping speed          | [arcmin^-2/hour/mJy^-2 ] | %6.1f  |'  % (s_map_spec) 
print '  ---------------------------------------------------------------'
print ' '

if verbose_flag == 1:

    print bcolors.WARNING + '  WARNING1: For proposed maps taking longer than 40 minutes, the total' 
    print bcolors.WARNING + '  integration time should be split into several maps.\n' + bcolors.ENDC
    print ' '
    print bcolors.WARNING + '  WARNING2: Note that the confusion noise is not included in this Time' 
    print bcolors.WARNING + '  Estimator and it should be discussed in the proposal whether appro-'
    print bcolors.WARNING + '  priate. The RMS noise does not decrease indefinitely with integra- '
    print bcolors.WARNING + '  tion time but stops at the confusion limit caused by unresolved sour-'
    print bcolors.WARNING + '  ces within the beam. The exact threshold at which the RMS noise    '
    print bcolors.WARNING + '  reaches the confusion limit will vary with wavelength, beam size,  '
    print bcolors.WARNING + '  and the type of astronomical source (Galactic or Extra-Galactic). For' 
    print bcolors.WARNING + '  instance, for the GOODS-N field (part of the Deep Field GT Proposal),' 
    print bcolors.WARNING + '  a preliminary 1 sigma confusion limit around ~ 0.090 mJy and ~ 0.056' 
    print bcolors.WARNING + '  mJy at 1.2 mm and 2.0 mm, respectively, has been estimated (A. Beelen,'
    print bcolors.WARNING + '  private communication).\n' + bcolors.ENDC

    print ' '
    print '             Time per OTF line for allowed scanning speeds   '
    print '        -----------------------------------------------------'
    print '        |          | OTF Scan Velocity  | TIME PER OTF LINE |'
    print '        -----------------------------------------------------'
    print '        |          |     [arcsec/s]     |    [minutes]      |'
    print '        -----------------------------------------------------'
    print '        | OTF-SLOW |       %6.1f       |     %6.2f        |' % (sampl_vel_slow,T_SLOW)
    print '        | OTF-MID  |       %6.1f       |     %6.2f        |' % (sampl_vel_mid,T_MID)
    print '        | OTF-FAST |       %6.1f       |     %6.2f        |' % (sampl_vel_fast,T_FAST)
    print '        -----------------------------------------------------'
    print ' '

###############################
#### CONVERSION TO MINUTES ####
###############################

if t_spec >= 1.:
   time_str_spec=str(round(t_spec,1))+' hours'

if t_spec < 1. and t_spec >= 1./60.:
   time_str_spec=str(round(t_spec*60,1))+' minutes'

if t_spec < 1./60. and t_spec >= 0.1/3600.:
   time_str_spec=str(int(t_spec*3600))+' seconds'

if t_spec < 0.1/3600.:
   time_str_spec='0.1 (seconds)'

slew_overhe = str(int((slewing_loss-1.0)*100.0))

##############################
#### OUTPUT FINAL RESULTS ####
##############################
    
#if dx < 6.5 or dy < 6.5: 
print bcolors.HEADER + bcolors.BOLD  + 'To optimize the correction of data instabilities, the size of scans along'
print bcolors.HEADER + bcolors.BOLD  + 'the scan direction should be at least:'  
print bcolors.HEADER + bcolors.BOLD  + 'NIKA2 FOV (6.5 arcminutes) + 2 * NIKA2 beam width (12 arcseconds at 1mm,' 
print bcolors.HEADER + bcolors.BOLD  + '18 arcseconds at 2mm) + source size above the noise + 2s * scan speed'
print bcolors.HEADER + bcolors.BOLD  +  bcolors.UNDERLINE + 'Please check your source size and your scanning speed to evaluate' 
print bcolors.HEADER + bcolors.BOLD  +  bcolors.UNDERLINE + 'the correct map size fitting your needs.' + bcolors.ENDC

print ''
print bcolors.HEADER + '***********************************************'
print '***  Total Integration Time => ' + time_str_spec + '   ***'
print '***********************************************\n'    
print ''

print bcolors.HEADER +  bcolors.UNDERLINE + 'Please include the following text into your proposal:\n' + bcolors.ENDC

print 'According to the published commissioning results of the NIKA2 instrument, the total observing time using the NIKA2 '+str(band)+' mm band to map a region of '+str(round(area_obs,1))+' [arcmin^2] to reach an rms of '+str(rms)+' [mJy/beam], assuming '+str(pwv)+' [mm] pwv, '+str(elevation)+' [deg] elevation, Filter = '+str(filtering)+', Overhead = '+str(foverhead)+', was estimated to be *'+time_str_spec+'*, using the time estimator '+version+'.\n'

