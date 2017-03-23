;; =========================NIKA & XPOL MAPS========================================
;; dir = '/Users/ritacco/Nika/Plots/CRAB_RADEC_LkgCorr'
;; map_1mm = mrdfits( dir+"/MAPS_1_Crab_v3.fits", 1, head_1mm)
;; map_2mm = mrdfits( dir+"/MAPS_2_Crab_v3.fits", 1, head_2mm)
dir = '/Users/ritacco/Nika/software/Processing/Labtools/AR/codes4thesis/CRAB/'
nika_ar  = mrdfits(dir+'ale2h_reso_xpol_Crab_2mm.fits', 1, head_2mm)
nika_stddev = mrdfits(dir+'ale2h_reso_xpol_Crab_2mm.fits', 2, head_2mm_stddev)

map_i_2mm = nika_ar[*,*,0]
map_q_2mm = nika_ar[*,*,1]
map_u_2mm = nika_ar[*,*,2]

map_i_stddev = nika_stddev[*,*,0]
map_q_stddev = nika_stddev[*,*,1]
map_u_stddev = nika_stddev[*,*,2]

map_i_3mm_ja = mrdfits('/Users/ritacco/Desktop/maps/crab_I_coadd28_reproj_nbl_ccentonly.fits', 0, h_ja)
map_q_3mm_ja = mrdfits('/Users/ritacco/Desktop/maps/crab_Q_coadd28_reproj_nbl_ccentonly.fits', 0, h_ja_q)
map_u_3mm_ja = mrdfits('/Users/ritacco/Desktop/maps/crab_U_coadd28_reproj_nbl_ccentonly.fits', 0, h_ja_u)

;; ========================= Adapting the astrometry ================================
extast,h_ja,astr
longobj = sxpar(h_ja,'crval1') 
latobj = sxpar(h_ja,'crval2') 
ra_xpol = SIXTY(longobj/15.0)
dec_xpol = SIXTY(latobj)
dec_xpol[2] = Float(Round(dec_xpol[2]*1000)/1000.)

nx_xpol = astr.naxis[0]
ny_xpol = astr.naxis[1]
x_xpol  = lindgen(nx_xpol)#replicate(1., ny_xpol)
y_xpol  = replicate(1,nx_xpol)#lindgen(ny_xpol)

astr.ctype[0] = "RA---TAN"
astr.ctype[1] = "DEC--TAN"
sxaddpar, h_ja, "naxis", 2
sxdelpar, h_ja, "naxis3"
sxaddpar, h_ja, "ctype1", "RA---TAN"
sxaddpar, h_ja, "ctype2", "DEC--TAN"

ad2xy, longobj,latobj, astr, x1, y1
xy2ad, x_xpol, y_xpol, astr, a, d
gcirc,2,longobj,latobj,a,d,dis
w_centre_xpol = where(dis/60. le 3.8,comp=w1, nw)

ra  = sxpar(h_ja, 'crval1')
dec = sxpar(h_ja, 'crval2')
nx  = sxpar(h_ja, 'naxis1')
ny  = sxpar(h_ja, 'naxis2')
reso= sxpar(h_ja, "CDELT2", /silent) * 3600.d0
naxis = [nx, ny]                                       ;Nb pixel along x and y
cd = [[1.0,-0.0],[0.0,1.0]]                            ;Rotation matrix but no rotation here
cdelt = [-1.0, 1.0] * reso/3600.0                      ;Pixel size (ra along -1)
crpix = [37.0,20.0]                                    ;Ref pixel (central pixel (always odd nb))
crval = [ra, dec]                                      ;ra dec of the ref pix
ctype = ['RA---TAN','DEC--TAN']                        ;Projection type
     
astrometry = {naxis:naxis, $
                   cd:cd, $
                   cdelt:cdelt, $
                   crpix:crpix, $
                   crval:crval, $
                   ctype:ctype,$
                   longpole:180.0, $
                   latpole:90.0, $
                   pv2:[0.0,0.0]}


mkhdr,  head_xpol, map_i_3mm_ja ;get header typique
putast, head_xpol, astrometry, equinox=2000, cd_type=0 ;astrometry in header
;; stop
extast,head_xpol,astr
longobj = sxpar(h_ja,'crval1') 
latobj  = sxpar(h_ja,'crval2') 
ra_xpol = SIXTY(longobj/15.0)
dec_xpol = SIXTY(latobj)
dec_xpol[2] = Float(Round(dec_xpol[2]*1000)/1000.)
coord_plot = [ten(ra_xpol[0],ra_xpol[1],ra_xpol[2])*15.0,$
              ten(dec_xpol[0],dec_xpol[1],dec_xpol[2])]    

;; ******************** Match Xpol and NIKA maps **************************
map_i = reform(nika_ar[*,*,0])
map_q = reform(nika_ar[*,*,1])
map_u = reform(nika_ar[*,*,2])
map_i_stddev = reform(nika_stddev[*,*,0])
map_q_stddev = reform(nika_stddev[*,*,1])
map_u_stddev = reform(nika_stddev[*,*,2])
omega_b_2mm = 2*!pi*[18.2*!fwhm2sigma]^2*corr
omega_pix = reso*reso

delta_x = -0.9
delta_y = 4.
corr = 1.28
map_i_nika = shift( map_i, delta_x, delta_y)
map_q_nika = shift( map_q, delta_x, delta_y)
map_u_nika = shift( map_u, delta_x, delta_y)
map_i_nika_stddev = shift( map_i_stddev, delta_x, delta_y)
map_q_nika_stddev = shift( map_q_stddev, delta_x, delta_y)
map_u_nika_stddev = shift( map_u_stddev, delta_x, delta_y)
col_cor = 1.05
q_nika_sigma2=map_q_nika_stddev^2*corr*(omega_pix/omega_b_2mm)^2*col_cor
u_nika_sigma2=map_u_nika_stddev^2*corr*(omega_pix/omega_b_2mm)^2*col_cor
i_nika_sigma2=map_i_nika_stddev^2*corr*(omega_pix/omega_b_2mm)^2*col_cor


ptg_offet = sqrt(delta_x^2+delta_y^2)*abs(sxpar(head_xpol, "cdelt1", /silent))*3600. ; arcsec
print, "Pointing offset: ", ptg_offet
map_ja   = map_i_3mm_ja*6. ;; Jy
;; map_ja[w1] = 0.d0
;; map_i_nika[w1] = 0.d0
cross_xpol_nika  = map_ja*map_i_nika

;; fov_plot = 500
;; overplot_radec_bar_map, map_i_nika, head_xpol, map_i_3mm_ja, head_xpol, fov_plot, reso, coord_plot,$
;;                         ;; pdf='/Volumes/STIPU/NIKA/Plots_tesi/crab_xpol_match2.pdf',$
;;                         conts1=[0.5, 1., 1.5, 2., 2.5], conts2=[0.2, 0.4, 0.8, 0.9], colconts2=0, $
;;                         range=[0.,2.8], colconts1=100, thickconts2=2,  /type, bg1=1e5, $
;;                         barcharthick=5, mapcharthick=5, barcharsize=1.5, mapcharsize=1. ,$
;;                         xtitle='R.A.', ytitle='Dec.'

print, total(cross_xpol_nika)

;; ==================== Covert Xpol maps in MJy/sr ===========================
ipol_xpol = sqrt(map_q_3mm_ja^2+map_u_3mm_ja^2)
lambda_micron = (300./90.)*10.^3
nk_convert_k2Mjy, lambda_micron, map_i_3mm_ja, sig_MJy_i_ja, /RJ
nk_convert_k2Mjy, lambda_micron, map_q_3mm_ja, sig_MJy_q_ja, /RJ
nk_convert_k2Mjy, lambda_micron, map_u_3mm_ja, sig_MJy_u_ja, /RJ
nk_convert_k2Mjy, lambda_micron, ipol_xpol, sig_MJy_ipol_ja, /RJ

int_rad = dindgen(360)/3600.
phi_i_ja     = nk_intgmap(sig_MJy_i_ja, 13.7*!arcsec2rad, int_rad*!dtor)*1e6
phi_q_ja     = nk_intgmap(sig_MJy_q_ja, 13.7*!arcsec2rad, int_rad*!dtor)*1e6
phi_u_ja     = nk_intgmap(sig_MJy_u_ja, 13.7*!arcsec2rad, int_rad*!dtor)*1e6
phi_ipol_ja  = nk_intgmap(sig_MJy_ipol_ja-0.5, 13.7*!arcsec2rad, int_rad*!dtor)*1e6

map_ja      = sig_MJy_i_ja*(13.7*!arcsec2rad)^2*1e6
map_ipol_ja = sig_MJy_ipol_ja*(13.7*!arcsec2rad)^2*1e6 ;; XPOL 27 arcsec FWHM 
;; ===================== Convert NIKA maps in Jy/sr ==========================
sxaddpar, head_2mm, "naxis", 2
sxdelpar, head_2mm, "naxis3"
sxaddpar, head_2mm, "ctype1", "RA---TAN"
sxaddpar, head_2mm, "ctype2", "DEC--TAN"
extast,head_2mm,astr
reso = astr.cdelt[1]*3600.0

thres= 3
nx   = sxpar(head_2mm, 'naxis1')
ny   = sxpar(head_2mm, 'naxis2')
xmap = reso*(replicate(1, ny) ## dindgen(nx)) - reso*(nx-1)/2.0
ymap = reso*(replicate(1, nx) #  dindgen(ny)) - reso*(ny-1)/2.0 
rmap = sqrt(xmap^2 + ymap^2)  
loc0 = where(rmap ge 120 and rmap ge 150 , nloc0)
omega_b_2mm = 2*!pi*(18.2*!fwhm2sigma)^2

zl2mm = mean(map_i_nika[loc0])
print, zl2mm
stop
corr=1.28
omega_pix = reso*reso
s = where( map_i_nika/map_i_nika_stddev gt thres and map_i_nika gt 0,comp=s1)
map_nika_baseline = map_i_nika; - zl2mm
;;i_2mm_nika = map_i_nika*omega_pix/(omega_b_2mm*corr)
i_2mm = map_nika_baseline*omega_pix/(omega_b_2mm*corr)
print, total(i_2mm)
ipol_2mm = sqrt(map_q_nika^2+map_u_nika^2 - map_q_nika_stddev^2 - map_u_nika_stddev^2 )
zl2mm_pol= mean(ipol_2mm[loc0])
ip_2mm   = ipol_2mm*omega_pix/(omega_b_2mm*corr)

header = head_xpol
fov_plot = 1000


;; define region of the crab
m_nika_i = i_2mm
m_xpol_i = map_ja
m_xpol_i[w1] = 0.d0
m_nika_i[w1] = 0.d0
print, total(map_ja)
print, total(i_2mm)

w_c = where(map_ja ge 0.4, comp=w2, nw)
m_x = map_ja
m_x[w2] = 0.d0
m_n = i_2mm
m_n[w2] = 0.d0


beta_nikaxpol = alog(total(m_x)/total(m_n))/alog(90.d0/150.d0)
print, beta_nikaxpol

w_c = where(map_ja ge 0.4, comp=w2, nw)
m_x = map_ja
m_x[w2] = 0.d0
m_n = i_2mm
m_x[w2] = 0.d0

beta_map = alog(m_x/m_n)/alog(90.d0/150.d0)
imview, beta_map, imr=[-0.4, 1]


;; overplot_radec_bar_map, map_ipol_ja, head_xpol, ip_2mm, header, fov_plot, reso, coord_plot,$
;;                         ;; pdf='/Volumes/STIPU/NIKA/Plots_tesi/crab_xpol_match2.pdf',$
;;                         conts1=[1e-5,1e5], conts2=[0.05, 0.1, 0.2, 0.28], colconts2=0, $
;;                         range=[0.,0.3], colconts1=100, thickconts2=2,  /type, bg1=1e5, $
;;                         barcharthick=5, mapcharthick=5, barcharsize=1.5, mapcharsize=1. ,$
;;                         xtitle='R.A.', ytitle='Dec.'

;; imview, map_ipol_ja, header=head_xpol, contour=ip_2mm, cont_header=header
;; ===================== Spectral index ======================================
;; COEFF = ROBUST_LINEFIT( alog(map_ja[s]), alog(i_2mm[s]), YFIT, SIG, COEF_SIG)  

;; plot,  map_ja[s], i_2mm[s], /xlog, /ylog, /nodata
;; oplot, map_ja[s], i_2mm[s], col=70, psym=8, symsize=0.2
;; oplot, map_ja[s], exp(yfit), col=230

;; m_xpol = map_ja
;; m_nika = i_2mm
;; m_xpol[s1] = !values.d_nan
;; m_nika[s1] = !values.d_nan
;; coeff0 = coeff[0] + coef_sig[1]*randomn(seed,100000)
;; beta_i = coeff[0]/alog(150./90.)
;; sigma_beta_i = stddev(coeff0)
;; print, 'beta fit'
;; print, beta_i, sigma_beta_i

;; map_ja[s1] = 0.d0
;; i_2mm[ s1] = 0.d0
;; spec_i   = [alog(map_ja/i_2mm)/alog(90./150.)]
;; spec_i[s1]=!values.d_nan
;; imview, spec_i, header=head_xpol, c_thick=5.,$
;;         c_colors=0, imr=[-0.05,1]
;; stop
;; set_plot, 'ps'
;; device,/color, bits_per_pixel=256, filename='/Users/ritacco/Documents/Plots/Crab/I_2mm_vs_3mm.ps'
;; plot,   map_ja[w_centre_xpol], i_2mm[w_centre_xpol], xtitle='90 GHz Stokes I [Jy/sr]', ytitle='150 GHz Stokes I [Jy/sr]',  thick=3, /xthick, /ythick, charsize=1.3, charthick=3., psym=8, symsize=0.2, /nodata, /xlog, /ylog
;; oplot,  map_ja[w_centre_xpol], i_2mm[w_centre_xpol], col=50, psym=8, symsize=0.2
;; oplot,  map_ja[w_centre_xpol], exp(yfit), col=230
;; legendastro, ['Best fit: '+string(coeff[1],format='(F5.2)')+' x !12x!3!U'+string(coeff[0],format='(F5.2)')], charsize=1.2, box=0, thick=6., charthick=3, /bottom
;; device,/close
;; set_plot, 'x'
;; pstopdf_crop, '/Users/ritacco/Documents/Plots/Crab/I_2mm_vs_3mm'
;; spawn, "\cp /Users/ritacco/Documents/Plots/Crab/I_2mm_vs_3mm.pdf /Users/ritacco/Nika/software/Processing/Papers_Proposal_Notes/Paper_crab/figures/."



levels = [0.5, 1., 1.5, 2., 2.5]
;; imview, spec_i,  header=head_xpol, contour=map_ja, imr=[-1,0.5],$
;;         cont_header=head_xpol, coltable=39, levels=levels, $
;;         unitsbar="Beta", xtitle='R. A.', ytitle='Dec.',$
;;         postscript='/Users/ritacco/Documents/Plots/Crab/crab_beta_i_xpol.ps'
;; pstopdf_crop, '/Users/ritacco/Documents/Plots/Crab/crab_beta_i_xpol'
;; spawn, "\cp /Users/ritacco/Documents/Plots/Crab/crab_beta_i_xpol.pdf /Users/ritacco/Nika/software/Processing/Papers_Proposal_Notes/Paper_crab/figures/."

;; ------------- Polarization spectral index map -----------------
COEFF = ROBUST_LINEFIT( alog(map_ipol_ja[s]), alog(ip_2mm[s]), YFIT, SIG, COEF_SIG)  
coeff0    = coeff[0] + coef_sig[1]*randomn(seed,100000)
beta_ipol = coeff[0]/alog(150./90.)
sigma_beta_ip = stddev(coeff0)
print, 'beta ipol fit'
print, beta_ipol, sigma_beta_ip
spec_ip    = [alog((map_ipol_ja)/(ip_2mm))/alog(90./150.)]
spec_ip[s1]=!values.d_nan

;; set_plot, 'ps'
;; device,/color, bits_per_pixel=256, filename='/Users/ritacco/Documents/Plots/Crab/Ipol_2mm_vs_3mm.ps'
;; plot,   map_ipol_ja[s], ip_2mm[s], xtitle='90 GHz Polarization intensity [Jy/sr]', $
;;         ytitle='150 GHz Polarization intensity [Jy/sr]',  thick=3, /xthick, /ythick, charsize=1.3, $
;;         charthick=3., psym=8, symsize=0.2, /xlog, /ylog, /nodata
;; oplot,  map_ipol_ja[s], ip_2mm[s], col=50, psym=8, symsize=0.2
;; oplot,  map_ipol_ja[s], exp(yfit), col=230
;; legendastro, ['Best fit: '+string(coeff[1],format='(F5.2)')+' x !12x!3!U'+string(coeff[0], format='(F5.2)')], charsize=1.2, box=0, thick=6., charthick=3, /bottom, /right
;; device,/close

;; set_plot, 'x'
;; pstopdf_crop, '/Users/ritacco/Documents/Plots/Crab/Ipol_2mm_vs_3mm'
;; spawn, "\cp /Users/ritacco/Documents/Plots/Crab/Ipol_2mm_vs_3mm.pdf /Users/ritacco/Nika/software/Processing/Papers_Proposal_Notes/Paper_crab/figures/."
;; levels = [0.5, 1., 1.5, 2., 2.5]
;; imview, spec_ip,  header=head_xpol, contour=map_ja, imr=[-2,1], $
;;         cont_header=head_xpol, coltable=39, levels=levels, $
;;         unitsbar="Beta", xtitle='R. A.', ytitle='Dec.',$
;;         postscript='/Users/ritacco/Documents/Plots/Crab/crab_beta_ip_xpol.ps'
;; pstopdf_crop, '/Users/ritacco/Documents/Plots/Crab/crab_beta_ip_xpol'
;; spawn, "\cp /Users/ritacco/Documents/Plots/Crab/crab_beta_ip_xpol.pdf /Users/ritacco/Nika/software/Processing/Papers_Proposal_Notes/Paper_crab/figures/."


;; ==========================================================================

;; **************************** SCUPOL **************************************

;; dir_scuba = '/Users/ritacco/SCUBAPOL/'
dir = '/Users/ritacco/Nika/software/Processing/Labtools/AR/codes4thesis/CRAB/'
;; m850 = mrdfits(dir_scuba+'scupollegacy_crab_cube.fits', 0, h_scuba)
scupol_map = mrdfits(dir+'ale2h_reso_xpol_Crab_850.fits', 1, h_scuba)
scupol_map_err = mrdfits(dir+'ale2h_reso_xpol_Crab_850.fits', 2, h_scuba)

sxaddpar, h_scuba, "naxis", 2
sxdelpar, h_scuba, "naxis3"
sxaddpar, h_scuba, "ctype1", "RA---TAN"
sxaddpar, h_scuba, "ctype2", "DEC--TAN"

longobj = sxpar(h_scuba,'crval1') 
latobj = sxpar(h_scuba,'crval2') 
ra_scuba = SIXTY(longobj/15.0)
dec_scuba = SIXTY(latobj)
dec_scuba[2] = Float(Round(dec_scuba[2]*1000)/1000.)
message, /info, 'IMBFITS pointing is :'
message, /info, 'R.A.: '+strtrim(ra_scuba[0],2)+' h '+strtrim(ra_scuba[1],2)+' min '+strtrim(ra_scuba[2],2)+' s'
message, /info, 'Dec.: '+strtrim(dec_scuba[0],2)+' deg '+strtrim(dec_scuba[1],2)+" arcmin "+strtrim(dec_scuba[2],2)+' arcsec'

coord_plot_scuba = [ten(ra_scuba[0],ra_scuba[1],ra_scuba[2])*15.0,$
                    ten(dec_scuba[0],dec_scuba[1],dec_scuba[2])]    

omega_scuba_850 = 2*!pi*(20.*!fwhm2sigma)^2

delta_x1 = 3
delta_y1 = 5.


mi = scupol_map[*,*,0]
mq = scupol_map[*,*,1]
mu = scupol_map[*,*,2]

mi_err = scupol_map_err[*,*,0]
mq_err = scupol_map_err[*,*,1]
mu_err = scupol_map_err[*,*,2]

map_scuba_i = shift( mi, delta_x1, delta_y1)
map_scuba_q = shift( mq, delta_x1, delta_y1)
map_scuba_u = shift( mu, delta_x1, delta_y1)

map_scuba_i_stddev = shift( mi_err, delta_x1, delta_y1)
map_scuba_q_stddev = shift( mq_err, delta_x1, delta_y1)
map_scuba_u_stddev = shift( mu_err, delta_x1, delta_y1)


;; fov_plot = 800
;; overplot_radec_bar_map, map_scuba, head_xpol, map_ja, head_xpol, fov_plot, reso, coord_plot,$
;;                         ;; pdf='/Volumes/STIPU/NIKA/Plots_tesi/crab_xpol_match2.pdf',$
;;                         conts1=[1e-5,1e5], conts2=[0.05, 0.1, 0.2, 0.9], colconts2=100, $
;;                         range=[0.,2.5], colconts1=100, thickconts2=2,  /type, bg1=1e5, $
;;                         barcharthick=5, mapcharthick=5, barcharsize=1.5, mapcharsize=1. ,$
;;                         xtitle='R.A.', ytitle='Dec.'

;; m850 = scupol_map[*,*,0]*omega_pix/omega_scuba_850
m850_i = map_scuba_i*omega_pix/omega_scuba_850
m850_q = map_scuba_q*omega_pix/omega_scuba_850
m850_u = map_scuba_u*omega_pix/omega_scuba_850
m850_i_stddev = map_scuba_i_stddev*omega_pix/omega_scuba_850
m850_q_stddev = map_scuba_q_stddev*omega_pix/omega_scuba_850
m850_u_stddev = map_scuba_u_stddev*omega_pix/omega_scuba_850

m850_i_sigma2 = total(map_scuba_i_stddev^2)*omega_pix/omega_scuba_850
m850_q_sigma2 = total(map_scuba_q_stddev^2)*omega_pix/omega_scuba_850
m850_u_sigma2 = total(map_scuba_u_stddev^2)*omega_pix/omega_scuba_850

iscuba_sigma2=total(m850_i_stddev^2)*(omega_pix/omega_scuba_850)^2
ipol_850 = sqrt(total(m850_q)^2+total(m850_u)^2 - total(m850_q_stddev^2)  - total(m850_u_stddev^2))
ipolscuba_err = sqrt(total(m850_q)^2*m850_q_sigma2+total(m850_u)^2*m850_u_sigma2)/(ipol_850)
print, 'Stokes I and error for SCUBA unit [Jy/sr]'
print, total(m850_i), iscuba_sigma2

print, 'IPOL and error for SCUBA unit [Jy/sr]'
print, ipol_850, ipolscuba_err
ip_850 = sqrt(m850_q^2+m850_u^2 - m850_q_stddev^2  - m850_u_stddev^2)

nu1 = !const.c/(850.*1e-6)/1e9
nu2 = !const.c/(!nika.lambda[1]*1e-3)/1e9

beta_im = alog(ip_850/ip_2mm)/alog(nu1/nu2)


stop
;; ipol_850 = sqrt(m850_q^2+m850_u^2)
beta = -0.294
ampl = 973+19

;; ******************convert maps to Jy and give the flux*************************
map_xpol_850micron   = sig_Mjy_i_ja*(nu1/90.)^beta
map_xpol_q_850micron = sig_Mjy_q_ja*(nu1/90.)^beta
map_xpol_u_850micron = sig_Mjy_u_ja*(nu1/90.)^beta
map_xpol_2mm         = sig_Mjy_i_ja*(nu2/90.)^beta
map_xpol_2           = map_xpol_2mm*omega_pix/omega_b_2mm

w = where(m850_i gt 0, comp=w1)
int_rad = dindgen(360)/3600.
phi_i_ja_850  = nk_intgmap(map_xpol_850micron[w], 13.7*!arcsec2rad, int_rad*!dtor)*1e6
;phi_i_scuba   = nk_intgmap(m850_i[w], reso, int_rad, var_map=m850_i_stddev^2, err=err_850)
phi_i_scuba   = nk_intgmap(m850_i[w], 13.7*!arcsec2rad, int_rad*!dtor)

map_2mm_scuba = m850_i[w]*(nu2/nu1)^beta
phi_i_ja_2mm  = nk_intgmap(map_xpol_2mm[w], 13.7*!arcsec2rad, int_rad*!dtor)*1e6
print, 'flux extrapolated at 2mm using the XPOL map'
print, phi_i_ja_2mm[n_elements(phi_i_ja_2mm)-1]
print, 'flux extrapolated at 2mm  using the scuba map'
print, total(map_2mm_scuba[w])
print, 'flux at 2mm measured by NIKA in the size SCUBA'
print, total(i_2mm[w])
print, 'flux measured by SCUBA'


beta_NIKA_scubapol_i = alog(total(m850_i[w])/total(i_2mm[w]))/alog(nu1/nu2)
beta_NIKA_scubapol_ip= alog(total(ip_850[w])/total(ip_2mm[w]))/alog(nu1/nu2)

beta_image = alog(total(ip_850)/total(ip_2mm))/alog(nu1/nu2)
beta_image[s1] = 0.d0

nika_intensity = total(i_2mm[w])
nika_err_i     = 12.95
nika_frequency = 150.d0

scuba_intensity= 98.3
scuba_err      = 0.1
scuba_frequency= 352.d0

i150  = nika_intensity+nika_err_i*randomn(seed,10000)
i850  = scuba_intensity+scuba_err*randomn(seed,10000)

betasim_i = alog(i850/i150)/alog(nu1/nu2)

print, 'beta flux'
print, mean(betasim_i), stddev(betasim_i)

;; error derivation for polarization intensity SED

nika_intensity_pol = total(ip_2mm[w])
nika_err_i     = 0.02
nika_frequency = 150.d0

scuba_intensity_pol= total(ip_850[w])
scuba_err      = 0.1
scuba_frequency= 352.d0

i150_ipol  = nika_intensity_pol+nika_err_i*randomn(seed,10000)
i850_ipol  = scuba_intensity_pol+scuba_err*randomn(seed,10000)

betasim_ipol = alog(i850_ipol/i150_ipol)/alog(nu1/nu2)

print, 'beta flux ipol'
print, mean(betasim_ipol), stddev(betasim_ipol)


print, 'spectral index considering SCUBA and NIKA data'
print, beta_NIKA_scubapol_i, stddev(betasim)
print, 'spectral index ipol SCUBA and NIKA data'
print, beta_NIKA_scubapol_ip, stddev(betasim_ipol)


stop

beta_NIKA_xpol = alog(total(map_ja[w])/total(i_2mm[w]))/alog(90./nu2)
print, 'spectral index considering XPOL and NIKA data'
print, beta_NIKA_xpol

col_cor=1.05
qscuba_nika = total(map_q_nika[w])*corr*omega_pix/omega_b_2mm*col_cor
uscuba_nika = total(map_u_nika[w])*corr*omega_pix/omega_b_2mm*col_cor
qscuba_nika_sigma2=total(map_q_nika_stddev[w]^2)*corr*(omega_pix/omega_b_2mm)^2*col_cor
uscuba_nika_sigma2=total(map_u_nika_stddev[w]^2)*corr*(omega_pix/omega_b_2mm)^2*col_cor
iscuba_nika = total(i_2mm[w])*1.05
iscuba_nika_sigma2=total(map_i_nika_stddev[w]^2)*corr*(omega_pix/omega_b_2mm)^2*col_cor

p = sqrt(qscuba_nika^2+uscuba_nika^2 - qscuba_nika_sigma2 - uscuba_nika_sigma2)/iscuba_nika
errp = sqrt(qscuba_nika^2*qscuba_nika_sigma2+uscuba_nika^2*uscuba_nika_sigma2 + p^4*iscuba_nika^2*iscuba_nika_sigma2)/(p*iscuba_nika^2)
print, 'pol_deg @ 1mm in scuba size'
print, p, errp

psi = 0.5*atan(uscuba_nika,qscuba_nika)*!radeg
sigma_psi = 0.5/(qscuba_nika^2 + uscuba_nika^2)*sqrt(qscuba_nika^2*qscuba_nika_sigma2^2 + uscuba_nika^2*uscuba_nika_sigma2^2)*!radeg
print, 'angle @ 1mm in scuba size'
print, psi+180, sigma_psi

ipolscuba_nika     = sqrt(qscuba_nika^2+uscuba_nika^2 - qscuba_nika_sigma2 - uscuba_nika_sigma2)
ipolscuba_nika_err = sqrt(qscuba_nika^2*qscuba_nika_sigma2+uscuba_nika^2*uscuba_nika_sigma2)/(ipolscuba_nika)
print, ' IPOL @ 1mm in scuba size'
print, ipolscuba_nika, ipolscuba_nika_err

print, 'Stokes I @ 1mm in scuba size'
print, iscuba_nika, sqrt(iscuba_nika_sigma2)




;;++++++++++++check nika and xpol maps in the scubapol size+

;; m_test_xpol  = map_ja
;; m_test_scuba = m850
;; m_test_nika  = i_2mm
;; m_test_xpol[w1]  = !values.d_nan
;; m_test_scuba[w1] = !values.d_nan
;; m_test_nika[w1 ] = !values.d_nan

;; wind,1
;; imview, m_test_xpol,  imr=[0,1.1],title = 'Map XPOL 90 GHz FoV of SCUBAPOL'
;; wind,2                
;; imview, m_test_scuba, imr=[0,1.1],title = 'Map SCUBAPOL 350 GHz'
;; wind,3                
;; imview, m_test_nika,  imr=[0,1.1],title = 'Map NIKA 150 GHz FoV SCUBAPOL'




stop

end



