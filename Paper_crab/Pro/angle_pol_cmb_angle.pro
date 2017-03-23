;; Calculating the Planck's pol angles in Ra,Dec coordinates
;; 30 GHz
glon = 184.543
glat = -5.782
I = 344.23 ;; Jy
sigmaI = 0.27
p = 7.10/100.
sigma = 0.33/100.
alpha = -89.26
sigma_alpha = 0.25 + 0.50
q = p*I*cos(2.d0*alpha*!dtor)
sigma_q= sqrt((I*sigma+p*sigmaI)^2*(cos(2.d0*alpha*!dtor))^2 + (p*I*2.d0*sin(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
u = p*I*sin(2.d0*alpha*!dtor)
sigma_u= sqrt((I*sigma+p*sigmaI)^2*(sin(2.d0*alpha*!dtor))^2 + (p*I*2.d0*cos(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
qu_gal2eq, glon, glat, q, u, q_out, u_out
alpha_deg = 0.5*atan(u_out,q_out)*!radeg 
sigma_alpha_deg = 0.5d0/(q_out^2+u_out^2)*sqrt( q_out^2*sigma_u^2 + u_out^2*sigma_q^2)*!radeg
print, 'polarization angle in radec coordinates at 30 GHz'
print, alpha_deg+180.d0
print, '+-', sigma_alpha_deg

;; 44 GHz
I = 292.68 ;; Jy
sigmaI = 0.23
p = 6.51/100.
sigma = 0.51/100.
alpha = -88.65
sigma_alpha = 0.79 + 0.50
q = p*I*cos(2.d0*alpha*!dtor)
u = p*I*sin(2.d0*alpha*!dtor)
qu_gal2eq, glon, glat, q, u, q_out, u_out
qu2alpha, q_out, u_out, alphadeg
sigma_q= sqrt((I*sigma+p*sigmaI)^2*(cos(2.d0*alpha*!dtor))^2 + (p*I*2.d0*sin(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_u= sqrt((I*sigma+p*sigmaI)^2*(sin(2.d0*alpha*!dtor))^2 + (p*I*2.d0*cos(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_alpha_deg = 0.5d0/(q_out^2+u_out^2)*sqrt( q_out^2*sigma_u^2 + u_out^2*sigma_q^2)*!radeg
print, 'polarization angle in radec coordinates at 44 GHz'
print, alphadeg+180.d0
print, '+-', sigma_alpha_deg


;; 70 GHz
I = 259.99 ;; Jy
sigmaI = 0.11
p = 7.9/100.
sigma = 0.32/100.
alpha = -87.49
sigma_alpha = 1.33 + 0.50
q = p*I*cos(2.d0*alpha*!dtor)
u = p*I*sin(2.d0*alpha*!dtor)
qu_gal2eq, glon, glat, q, u, q_out, u_out
qu2alpha, q_out, u_out, alphadeg
sigma_q= sqrt((I*sigma+p*sigmaI)^2*(cos(2.d0*alpha*!dtor))^2 + (p*I*2.d0*sin(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_u= sqrt((I*sigma+p*sigmaI)^2*(sin(2.d0*alpha*!dtor))^2 + (p*I*2.d0*cos(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_alpha_deg = 0.5d0/(q_out^2+u_out^2)*sqrt( q_out^2*sigma_u^2 + u_out^2*sigma_q^2)*!radeg
print, 'polarization angle in radec coordinates at 70 GHz'
print, alphadeg+180.d0
print, '+-', sigma_alpha_deg



;; 100 GHz
I = 215.16 ;; Jy
sigmaI = 0.06
p = 7.22/100.
sigma = 0.06/100.
alpha = -87.59
sigma_alpha = 0.26 + 0.62
q = p*I*cos(2.d0*alpha*!dtor)
u = p*I*sin(2.d0*alpha*!dtor)
qu_gal2eq, glon, glat, q, u, q_out, u_out
qu2alpha, q_out, u_out, alphadeg
sigma_q= sqrt((I*sigma+p*sigmaI)^2*(cos(2.d0*alpha*!dtor))^2 + (p*I*2.d0*sin(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_u= sqrt((I*sigma+p*sigmaI)^2*(sin(2.d0*alpha*!dtor))^2 + (p*I*2.d0*cos(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_alpha_deg = 0.5d0/(q_out^2+u_out^2)*sqrt( q_out^2*sigma_u^2 + u_out^2*sigma_q^2)*!radeg
print, 'polarization angle in radec coordinates at 100 GHz'
print, alphadeg+180.d0
print, '+-', sigma_alpha_deg


;; 143 GHz
I = 167.10 ;; Jy
sigmaI = 0.04
p = 7.19/100.
sigma = 0.05/100.
alpha = -87.03
sigma_alpha = 0.35 + 0.62
q = p*I*cos(2.d0*alpha*!dtor)
u = p*I*sin(2.d0*alpha*!dtor)
qu_gal2eq, glon, glat, q, u, q_out, u_out
qu2alpha, q_out, u_out, alphadeg
sigma_q= sqrt((I*sigma+p*sigmaI)^2*(cos(2.d0*alpha*!dtor))^2 + (p*I*2.d0*sin(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_u= sqrt((I*sigma+p*sigmaI)^2*(sin(2.d0*alpha*!dtor))^2 + (p*I*2.d0*cos(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_alpha_deg = 0.5d0/(q_out^2+u_out^2)*sqrt( q_out^2*sigma_u^2 + u_out^2*sigma_q^2)*!radeg
print, 'polarization angle in radec coordinates at 143 GHz'
print, alphadeg+180.d0
print, '+-', sigma_alpha_deg


;; 217 GHz
I = 124.21 ;; Jy
sigmaI = 0.04
p = 8.12/100.
sigma = 0.06/100.
alpha = -88.84
sigma_alpha = 0.55 + 0.62
q = p*I*cos(2.d0*alpha*!dtor)
u = p*I*sin(2.d0*alpha*!dtor)
qu_gal2eq, glon, glat, q, u, q_out, u_out
qu2alpha, q_out, u_out, alphadeg
sigma_q= sqrt((I*sigma+p*sigmaI)^2*(cos(2.d0*alpha*!dtor))^2 + (p*I*2.d0*sin(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_u= sqrt((I*sigma+p*sigmaI)^2*(sin(2.d0*alpha*!dtor))^2 + (p*I*2.d0*cos(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_alpha_deg = 0.5d0/(q_out^2+u_out^2)*sqrt( q_out^2*sigma_u^2 + u_out^2*sigma_q^2)*!radeg
print, 'polarization angle in radec coordinates at 217 GHz'
print, alphadeg+180.d0
print, '+-', sigma_alpha_deg


;; 353 GHz
I = 82.17 ;; Jy
sigmaI = 0.67
p = 12.02/100.
sigma = 0.23/100.
alpha = -85.16
sigma_alpha = 1.93 + 0.62
q = p*I*cos(2.d0*alpha*!dtor)
u = p*I*sin(2.d0*alpha*!dtor)
qu_gal2eq, glon, glat, q, u, q_out, u_out
qu2alpha, q_out, u_out, alphadeg
sigma_q= sqrt((I*sigma+p*sigmaI)^2*(cos(2.d0*alpha*!dtor))^2 + (p*I*2.d0*sin(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_u= sqrt((I*sigma+p*sigmaI)^2*(sin(2.d0*alpha*!dtor))^2 + (p*I*2.d0*cos(2.d0*alpha*!dtor)*sigma_alpha*!dtor)^2)
sigma_alpha_deg = 0.5d0/(q_out^2+u_out^2)*sqrt( q_out^2*sigma_u^2 + u_out^2*sigma_q^2)*!radeg
print, 'polarization angle in radec coordinates at 353 GHz'
print, alphadeg+180.d0
print, '+-', sigma_alpha_deg



;;********************** check NIKA values in galactic coordinates ***************************
dir = '/Volumes/STIPU/NIKA/Plots/CRAB_RADEC_LkgCorr'
map_2mm = mrdfits( dir+"/MAPS_2_Crab_v3.fits", 1, head_2mm)
map_2mm_stddev = mrdfits( dir+"/MAPS_2_Crab_v3.fits", 2, header)
map_i_2mm = map_2mm[*,*,0]
map_q_2mm = map_2mm[*,*,1]
map_u_2mm = map_2mm[*,*,2]
map_i_2mm_stddev = map_2mm_stddev[*,*,0]
map_q_2mm_stddev = map_2mm_stddev[*,*,1]
map_u_2mm_stddev = map_2mm_stddev[*,*,2]
extast,head_2mm,astr

longobj = sxpar(head_2mm,'crval1') 
latobj = sxpar(head_2mm,'crval2') 
ra = SIXTY(longobj/15.0)
dec = SIXTY(latobj)
dec[2] = Float(Round(dec[2]*1000)/1000.)

nx = astr.naxis[0]
ny = astr.naxis[1]
x = lindgen(nx)#replicate(1., ny)
y = replicate(1,nx)#lindgen(ny)
reso = astr.cdelt[1]*3600.0
omega_pix = reso*reso
corr = 1.28
omega_b = 2*!pi*[18.2*!fwhm2sigma]^2*corr
col_cor = 1.05

astr.ctype[0] = "RA---TAN"
astr.ctype[1] = "DEC--TAN"
ad2xy, longobj,latobj, astr, x1, y1
xy2ad, x, y, astr, a, d
gcirc,2,longobj,latobj,a,d,dis
qu_eq2gal,a,d,map_q_2mm, map_u_2mm, Qgal_2mm,Ugal_2mm
qu_eq2gal,a,d,map_q_2mm_stddev, map_u_2mm_stddev, Qgal_2mm_stddev,Ugal_2mm_stddev

w = where(dis/60. le 5 and abs(map_q_2mm/map_q_2mm_stddev) gt 3 and abs(map_u_2mm/map_u_2mm_stddev) gt 3, comp=w1 , nw)
tot_i_2mm = total(map_i_2mm[w])*(omega_pix/omega_b )*col_cor 
tot_q_gal_2mm = total(Qgal_2mm[w])*col_cor*(omega_pix/omega_b)
tot_u_gal_2mm = total(Ugal_2mm[w])*col_cor*(omega_pix/omega_b)
tot_i_2mm_sigma2 = total((map_i_2mm_stddev[w])^2)*col_cor*(omega_pix/omega_b)^2
tot_q_gal_2mm_sigma2 = total((Qgal_2mm_stddev[w])^2)*col_cor*(omega_pix/omega_b)^2
tot_u_gal_2mm_sigma2 = total((Ugal_2mm_stddev[w])^2)*col_cor*(omega_pix/omega_b)^2
sigma_psi = 0.5/(tot_q_gal_2mm^2 + tot_u_gal_2mm^2)*sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2 + tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2)*!radeg
print, 'angle in galactic coordinates in 5 arcminutes and > 3 sigma'
psi = 0.5*atan(-tot_u_gal_2mm, tot_q_gal_2mm)*!radeg
print,  psi
print, sigma_psi
print, 'intensity', tot_i_2mm, sqrt(tot_i_2mm_sigma2)

p2mm = sqrt(tot_q_gal_2mm^2+tot_u_gal_2mm^2 - tot_q_gal_2mm_sigma2 - tot_q_gal_2mm_sigma2)/tot_i_2mm
errp2mm = sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2+tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2 + p2mm^4*tot_i_2mm^2*tot_i_2mm_sigma2)/(p2mm*tot_i_2mm^2)
print, 'pol deg', p2mm, errp2mm
ip   = sqrt(tot_q_gal_2mm^2 + tot_u_gal_2mm^2 - tot_q_gal_2mm_sigma2 - tot_u_gal_2mm_sigma2)
ip_err = sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2+tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2)/(ip)
print, 'pol intensity'
print, ip, ip_err

s = size(map_i_2mm)
nx = s[1]
ny = s[2]
p_2mm_sim   = dblarr(1000)
psi_2mm_sim = dblarr(1000)


for i=0, 999 do begin
   q_sim = p2mm[0]*map_i_2mm*cos(2.d0*(-87)*!dtor)+randomn(seed, nx, ny)*Qgal_2mm_stddev
   u_sim = p2mm[0]*map_i_2mm*sin(2.d0*(-87)*!dtor)+randomn(seed, nx, ny)*Ugal_2mm_stddev
   q5arc = total(q_sim[w])*omega_pix/omega_b*col_cor
   u5arc = total(u_sim[w])*omega_pix/omega_b*col_cor
   q5arc_sigma2=total(Qgal_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   u5arc_sigma2=total(Ugal_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   i5arc = total(map_i_2mm[w])*omega_pix/omega_b*col_cor
   i5arc_sigma2=total(map_i_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   p_2mm_sim[i]   = sqrt(q5arc^2+u5arc^2 - q5arc_sigma2 - u5arc_sigma2)/i5arc  
   psi_2mm_sim[i] = 0.5*atan(u5arc,q5arc)*!radeg
endfor
print, 'Polarization results obtained by using a simulation of the noise in 5 arcmin and P>3sigma'
print, 'angle'
print, mean(psi_2mm_sim), stddev(psi_2mm_sim)
print, 'polarization degree'
print, mean(p_2mm_sim)*100, stddev(p_2mm_sim)*100
stop

w = where(dis/60. le 5 and map_i_2mm_stddev gt 0 and map_q_2mm_stddev gt 0 and map_u_2mm_stddev gt 0, comp=w1)
tot_i_2mm = total(map_i_2mm[w])*(omega_pix/omega_b)*col_cor 
tot_q_gal_2mm = total(Qgal_2mm[w])*col_cor*(omega_pix/omega_b)
tot_u_gal_2mm = total(Ugal_2mm[w])*col_cor*(omega_pix/omega_b)
tot_i_2mm_sigma2 = total((map_i_2mm_stddev[w])^2)*col_cor*(omega_pix/omega_b)^2
tot_q_gal_2mm_sigma2 = total((Qgal_2mm_stddev[w])^2)*col_cor*(omega_pix/omega_b)^2
tot_u_gal_2mm_sigma2 = total((Ugal_2mm_stddev[w])^2)*col_cor*(omega_pix/omega_b)^2
sigma_psi = 0.5/(tot_q_gal_2mm^2 + tot_u_gal_2mm^2)*sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2 + tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2)*!radeg
print, 'angle in galactic coordinates in 5 arcminutes'
stop
print,  0.5*atan(-tot_u_gal_2mm,tot_q_gal_2mm)*!radeg
print, sigma_psi
print, 'intensity', tot_i_2mm, sqrt(tot_i_2mm_sigma2)
p2mm = sqrt(tot_q_gal_2mm^2+tot_u_gal_2mm^2 - tot_q_gal_2mm_sigma2 - tot_q_gal_2mm_sigma2)/tot_i_2mm
errp2mm = sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2+tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2 + p2mm^4*tot_i_2mm^2*tot_i_2mm_sigma2)/(p2mm*tot_i_2mm^2)
print, 'pol deg', p2mm, errp2mm
ip   = sqrt(tot_q_gal_2mm^2 + tot_u_gal_2mm^2 - tot_q_gal_2mm_sigma2 - tot_u_gal_2mm_sigma2)
ip_err = sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2+tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2)/(ip)
print, 'pol intensity'
print, ip, ip_err

for i=0, 999 do begin
   q_sim = p2mm[0]*map_i_2mm*cos(2.d0*(-87)*!dtor)+randomn(seed, nx, ny)*Qgal_2mm_stddev
   u_sim = p2mm[0]*map_i_2mm*sin(2.d0*(-87)*!dtor)+randomn(seed, nx, ny)*Ugal_2mm_stddev
   q5arc = total(q_sim[w])*omega_pix/omega_b*col_cor
   u5arc = total(u_sim[w])*omega_pix/omega_b*col_cor
   q5arc_sigma2=total(Qgal_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   u5arc_sigma2=total(Ugal_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   i5arc = total(map_i_2mm[w])*omega_pix/omega_b*col_cor
   i5arc_sigma2=total(map_i_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   p_2mm_sim[i]   = sqrt(q5arc^2+u5arc^2 - q5arc_sigma2 - u5arc_sigma2)/i5arc  
   psi_2mm_sim[i] = 0.5*atan(u5arc,q5arc)*!radeg
endfor
print, 'Polarization results obtained by using a simulation of the noise in 5 arcmin'
print, 'angle'
print, mean(psi_2mm_sim), stddev(psi_2mm_sim)
print, 'polarization degree'
print, mean(p_2mm_sim)*100, stddev(p_2mm_sim)*100
stop


w = where(dis/60. le 7 and map_i_2mm_stddev gt 0 and map_q_2mm_stddev gt 0 and map_u_2mm_stddev gt 0 and map_i_2mm gt 0)
tot_i_2mm = total(map_i_2mm[w])*col_cor*(omega_pix/omega_b)
tot_q_gal_2mm = total(Qgal_2mm[w])*col_cor*(omega_pix/omega_b)
tot_u_gal_2mm = total(Ugal_2mm[w])*col_cor*(omega_pix/omega_b)
tot_i_2mm_sigma2 = total((map_i_2mm_stddev[w])^2)*col_cor*(omega_pix/omega_b)^2
tot_q_gal_2mm_sigma2 = total((Qgal_2mm_stddev[w])^2)*col_cor*(omega_pix/omega_b)^2
tot_u_gal_2mm_sigma2 = total((Ugal_2mm_stddev[w])^2)*col_cor*(omega_pix/omega_b)^2
ip   = sqrt(tot_q_gal_2mm^2 + tot_u_gal_2mm^2 - tot_q_gal_2mm_sigma2 - tot_u_gal_2mm_sigma2)
ip_err = sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2+tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2)/(ip)
print, 'angle in galactic coordinates in 7 arcminutes'
print,  0.5*atan(-tot_u_gal_2mm,tot_q_gal_2mm)*!radeg

sigma_psi = 0.5/(tot_q_gal_2mm^2 + tot_u_gal_2mm^2)*sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2 + tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2)*!radeg
print, sigma_psi
print, 'intensity', tot_i_2mm, sqrt(tot_i_2mm_sigma2)
p2mm = sqrt(tot_q_gal_2mm^2+tot_u_gal_2mm^2 - tot_q_gal_2mm_sigma2 - tot_q_gal_2mm_sigma2)/tot_i_2mm
errp2mm = sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2+tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2 + p2mm^4*tot_i_2mm^2*tot_i_2mm_sigma2)/(p2mm*tot_i_2mm^2)
print, 'pol deg', p2mm, errp2mm
print, 'pol intensity'
print, ip, ip_err

for i=0, 999 do begin
   q_sim = p2mm[0]*map_i_2mm*cos(2.d0*(-87)*!dtor)+randomn(seed, nx, ny)*Qgal_2mm_stddev
   u_sim = p2mm[0]*map_i_2mm*sin(2.d0*(-87)*!dtor)+randomn(seed, nx, ny)*Ugal_2mm_stddev
   q5arc = total(q_sim[w])*omega_pix/omega_b*col_cor
   u5arc = total(u_sim[w])*omega_pix/omega_b*col_cor
   q5arc_sigma2=total(Qgal_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   u5arc_sigma2=total(Ugal_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   i5arc = total(map_i_2mm[w])*omega_pix/omega_b*col_cor
   i5arc_sigma2=total(map_i_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   p_2mm_sim[i]   = sqrt(q5arc^2+u5arc^2 - q5arc_sigma2 - u5arc_sigma2)/i5arc  
   psi_2mm_sim[i] = 0.5*atan(u5arc,q5arc)*!radeg
endfor
print, 'Polarization results obtained by using a simulation of the noise in 7 arcmin'
print, 'angle'
print, mean(psi_2mm_sim), stddev(psi_2mm_sim)
print, 'polarization degree'
print, mean(p_2mm_sim)*100, stddev(p_2mm_sim)*100

stop

ipol   = sqrt(Qgal_2mm^2+Ugal_2mm^2 -Qgal_2mm_stddev^2  - Ugal_2mm_stddev^2)

w = where(dis/60. le 1.4 and ipol gt 0.2, comp=w1 , nw)
m = ipol
m[w1] = 0.d0
imview, m

tot_i_2mm = total(map_i_2mm[w])*(omega_pix/omega_b)
tot_q_gal_2mm = total(Qgal_2mm[w])*(omega_pix/omega_b)
tot_u_gal_2mm = total(ugal_2mm[w])*(omega_pix/omega_b)
tot_q_gal_2mm_sigma2 = total((Qgal_2mm_stddev[w])^2)*(omega_pix/omega_b)^2
tot_u_gal_2mm_sigma2 = total((Ugal_2mm_stddev[w])^2)*(omega_pix/omega_b)^2
ip   = sqrt(tot_q_gal_2mm^2 + tot_u_gal_2mm^2 - tot_q_gal_2mm_sigma2 - tot_u_gal_2mm_sigma2)
ip_err = sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2+tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2)/(ip)
print, 'angle in galactic coordinates in 1.4 arcminutes and ipol > 0.2 Jy'
print,  0.5*atan(-tot_u_gal_2mm,tot_q_gal_2mm)*!radeg
sigma_psi = 0.5/(tot_q_gal_2mm^2 + tot_u_gal_2mm^2)*sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2 + tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2)*!radeg
print, sigma_psi
print, 'intensity', tot_i_2mm, sqrt(tot_i_2mm_sigma2)
p2mm = sqrt(tot_q_gal_2mm^2+tot_u_gal_2mm^2 - tot_q_gal_2mm_sigma2 - tot_q_gal_2mm_sigma2)/tot_i_2mm
errp2mm = sqrt(tot_q_gal_2mm^2*tot_q_gal_2mm_sigma2+tot_u_gal_2mm^2*tot_u_gal_2mm_sigma2 + p2mm^4*tot_i_2mm^2*tot_i_2mm_sigma2)/(p2mm*tot_i_2mm^2)
print, 'pol deg', p2mm, errp2mm
print, 'pol intensity'
print, ip, ip_err

for i=0, 999 do begin
   q_sim = p2mm[0]*map_i_2mm*cos(2.d0*(-87)*!dtor)+randomn(seed, nx, ny)*Qgal_2mm_stddev
   u_sim = p2mm[0]*map_i_2mm*sin(2.d0*(-87)*!dtor)+randomn(seed, nx, ny)*Ugal_2mm_stddev
   q5arc = total(q_sim[w])*omega_pix/omega_b*col_cor
   u5arc = total(u_sim[w])*omega_pix/omega_b*col_cor
   q5arc_sigma2=total(Qgal_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   u5arc_sigma2=total(Ugal_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   i5arc = total(map_i_2mm[w])*omega_pix/omega_b*col_cor
   i5arc_sigma2=total(map_i_2mm_stddev[w]^2)*(omega_pix/omega_b)^2*col_cor
   p_2mm_sim[i]   = sqrt(q5arc^2+u5arc^2 - q5arc_sigma2 - u5arc_sigma2)/i5arc  
   psi_2mm_sim[i] = 0.5*atan(u5arc,q5arc)*!radeg
endfor
print, 'Polarization results obtained by using a simulation of the noise in 1.4 arcmin'
print, 'angle'
print, mean(psi_2mm_sim), stddev(psi_2mm_sim)
print, 'polarization degree'
print, mean(p_2mm_sim)*100, stddev(p_2mm_sim)*100

stop

I = 84.6
p = 12.4/100.
alpha = -87.15
glon = 184.543
glat = -5.782
q = p*I*cos(2.d0*alpha*!dtor)
u = p*I*sin(2.d0*alpha*!dtor)
qu_gal2eq, glon, glat, q, u, q_out, u_out
;; qu_gal2eq, glon, glat, total(Qgal_2mm[w]), total(Ugal_2mm[w]), q_out, u_out
qu2alpha, q_out, u_out, alphadeg
print, alphadeg+180.d0
;;********************** check NIKA values in galactic coordinates ***************************

;; WMAP
q_wmap = -27.13
u_wmap = 1.4
alpha  = -88.5

qu_gal2eq, glon, glat, q_wmap, u_wmap, q_o, u_o
qu2alpha, q_o, u_o, alphadeg
print, alphadeg+180.d0



;; ++++++++++++++++++++++++++ SCUPOL results ++++++++++++++++++++++++++++++++++
scupol_map = mrdfits('/Users/ritacco/SCUBAPOL/scupollegacy_crab_cube.fits', 0, head_scuba)
scupol_map_err = mrdfits('/Users/ritacco/SCUBAPOL/scupollegacy_crab_cube.fits', 1, head_scuba) ;; variance

mi = scupol_map[*,*,0]*455.d0
mq = scupol_map[*,*,1]*455.d0
mu = scupol_map[*,*,2]*455.d0

mi_err = scupol_map_err[*,*,0]*455.d0
mq_err = scupol_map_err[*,*,1]*455.d0
mu_err = scupol_map_err[*,*,2]*455.d0
extast,head_scuba,astr
longobj = sxpar(head_scuba,'crval1') 
latobj = sxpar(head_scuba,'crval2') 
nx = astr.naxis[0]
ny = astr.naxis[1]
x = lindgen(nx)#replicate(1., ny)
y = replicate(1,nx)#lindgen(ny)
omega_pix = reso*reso
astr.ctype[0] = "RA---TAN"
astr.ctype[1] = "DEC--TAN"
ad2xy, longobj,latobj, astr, x1, y1
xy2ad, x, y, astr, a, d
gcirc,2,longobj,latobj,a,d,dis
qu_eq2gal,a,d,mq, mu, Qgal_scuba,Ugal_scuba
qu_eq2gal,a,d,mq_err, mu_err, Qgal_scuba_stddev,Ugal_scuba_stddev
omega_scuba_850 = 2*!pi*(20.*!fwhm2sigma)^2


ipol   = sqrt(Qgal_scuba^2+Ugal_scuba^2 -Qgal_scuba_stddev^2  - Ugal_scuba_stddev^2)
w = where(dis/60. le 1.4 and ipol gt 0.01, comp=w1 , nw)
w = where(ipol gt 0.2 and dis/60. le 1.4, comp=w1 , nw)
t = where(mi gt 0. and abs(mq) gt 0. and abs(mu) gt 0., comp=t1)
m = ipol
i = mi
m[w1] = 0.d0
i[w1] = 0.d0
imview, m
imview, i
tot_i_scuba = total(mi[w])*(omega_pix/omega_scuba_850)
tot_q_gal_scuba = total(Qgal_scuba[w])*omega_pix/omega_scuba_850
tot_u_gal_scuba = total(Ugal_scuba[w])*omega_pix/omega_scuba_850
tot_i_scuba_sigma2 = total((mi_err[w]))*(omega_pix/omega_scuba_850)^2
tot_q_gal_scuba_sigma2 = total((Qgal_scuba_stddev[w]))*(omega_pix/omega_scuba_850)^2
tot_u_gal_scuba_sigma2 = total((Ugal_scuba_stddev[w]))*(omega_pix/omega_scuba_850)^2
sigma_psi = 0.5/(tot_q_gal_scuba^2 + tot_u_gal_scuba^2)*sqrt(tot_q_gal_scuba^2*tot_q_gal_scuba_sigma2 + tot_u_gal_scuba^2*tot_u_gal_scuba_sigma2)*!radeg


print, 'angle in galactic coordinates in 1.4 arcminutes and > 3 sigma SCUBA'
print,  0.5*atan(-total(Ugal_scuba[w])*omega_pix/omega_scuba_850,total(Qgal_scuba[w])*omega_pix/omega_scuba_850)*!radeg
;; print,  0.5*atan(total(mu[w])*omega_pix/omega_scuba_850,total(mq[w])*omega_pix/omega_scuba_850)*!radeg+180.
print, sigma_psi
print, 'intensity', tot_i_scuba, sqrt(tot_i_scuba_sigma2)
pscuba = sqrt(tot_q_gal_scuba^2+tot_u_gal_scuba^2 - tot_q_gal_scuba_sigma2 - tot_u_gal_scuba_sigma2)/tot_i_scuba
errpscuba = sqrt(tot_q_gal_scuba^2*tot_q_gal_scuba_sigma2+tot_u_gal_scuba^2*tot_u_gal_scuba_sigma2 + pscuba^4*tot_i_scuba^2*tot_i_scuba_sigma2)/(pscuba*tot_i_scuba^2)
print, 'pol deg', pscuba, errpscuba
ip   = sqrt(tot_q_gal_scuba^2 + tot_u_gal_scuba^2 - tot_q_gal_scuba_sigma2 - tot_u_gal_scuba_sigma2)
ip_stddev = sqrt(tot_q_gal_scuba^2*tot_q_gal_scuba_sigma2+tot_u_gal_scuba^2*tot_u_gal_scuba_sigma2)/(ip)
print, 'pol intensity'
print, ip, ip_stddev


;; XPOL

dir_xpol = '/Users/ritacco/Desktop/maps/'

map_ja_i = mrdfits(dir_xpol+'crab_I_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
map_ja_q = mrdfits(dir_xpol+'crab_Q_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
map_ja_u = mrdfits(dir_xpol+'crab_U_coadd28_reproj_nbl_ccentonly.fits',0,head_ja)
extast,head_ja,astr
astr.ctype[0] = "RA---TAN"
astr.ctype[1] = "DEC--TAN"

reso = astr.cdelt[1]*3600.0

map_ja_i = map_ja_i*6.; Jy factor conversion Aumont et al.
map_ja_q = map_ja_q*6.; Jy factor conversion Aumont et al.
map_ja_u = map_ja_u*6.; Jy factor conversion Aumont et al.

longobj = sxpar(head_ja,'crval1') 
latobj = sxpar(head_ja,'crval2') 
nx = astr.naxis[0]
ny = astr.naxis[1]
x = lindgen(nx)#replicate(1., ny)
y = replicate(1,nx)#lindgen(ny)
omega_pix = reso*reso
ad2xy, longobj,latobj, astr, x1, y1
xy2ad, x, y, astr, a, d
gcirc,2,longobj,latobj,a,d,dis
qu_eq2gal,a,d,map_ja_q, map_ja_u, Qgal_ja,Ugal_ja


ipol   = sqrt(Qgal_ja^2+Ugal_ja^2)
w = where(dis/60. le 5 and map_ja_i gt 0 , comp=w1 , nw)
tot_i_ja = total(map_ja_i[w])
tot_q_gal_ja = total(Qgal_ja[w])
tot_u_gal_ja = total(Ugal_ja[w])
print,  0.5*atan(total(Ugal_ja[w]),total(Qgal_ja[w]))*!radeg
print, sqrt(total(Ugal_ja[w])^2+total(Qgal_ja[w])^2)/tot_i_ja

























end



