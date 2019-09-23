
;; Simple component separation to assess impact of non linearity

;; Assumed non-linearity parameter
epsilon = 1d-4 ;0.d0 ;1d-3 ;0.5d-1 ;1d-6 ;0.d0

;; Foreground mask parameter
latitude_cut = 30.

;; Outplot format
ps  = 0
png = 0
noplot = 0

;; Map size
map_xsize  = 20. ; 10.
map_ysize = map_xsize

;; Power of two to save time in FFT's
r = 9.
nx = 2L^r
ny = 2L^r
res_arcmin = map_xsize*60./nx ; arcmin

;; Foreground spectral parameters
beta_dust = 1.5
beta_sync = -3

;; Experiment bands
;; nu = [70.d0, 100.d0, 143.d0, 217.d0, 353.d0];, 545.d0]
nu = [70.d0, 143.d0, 353.d0]

;; Take a n-detector system with perfect angular coverage per band
alpha_det_deg = [0.d0, 60.d0, 120.d0]

;; MC realizations
nmc = 10 ; 50 ;30

;;============================================
n_nu  = n_elements(nu)
alpha = alpha_det_deg*!dtor
ndet  = n_elements(alpha)

;; Quick check
readcol, '$SK_DIR/Cl/cmb_totCls_r0.001.dat', l, clt, cle, clb, clte, $
         format='D,D,D,D', comment='#'
;; prefactor included in the CAMB's output files
fl = l*(l+1)/(2*!dpi)
clt  /= fl
cle  /= fl
clb_in = clb/fl
clte /= fl
;; No B in the simulations
clb = cle*0.d0
cls2mapsiqu, l, clt, cle, clb, clte, nx, res_arcmin/60., $
             cmb_t, cmb_q, cmb_u, cu_t, cu_e, cu_b, cu_te


;; Get foreground maps and derive fit power spectra
nside = 256
if defined(i_dust_ref) eq 0 then begin
   get_planck_maps, nside, $
                    i_cmb, q_cmb, u_cmb, $
                    i_dust_ref, q_dust_ref, u_dust_ref, $
                    i_sync_ref, q_sync_ref, u_sync_ref, /noplot
endif

npix = nside2npix(nside)
cover = dblarr(npix)+1.d0

;; Mask out the galactic plane (OK too)
if latitude_cut gt 0 then begin
   cover        = dblarr(npix)
   ipring       = lindgen(npix)
   pix2ang_ring, nside, ipring, theta, phi
   latitude = 90.-theta*!radeg
   w = where( abs(latitude) ge latitude_cut, nw)
   cover[w] = 1.d0
endif

dust_p_nu0 = 353.d0
dust_i_nu0 = 545.d0
sync_p_nu0 = 30.d0
sync_i_nu0 = 0.408d0

;; check at 95GHz to compare to Planck 2018, XI
nu0 = 95.

i_sync_ref *= (nu0/sync_i_nu0)^(-3)
q_sync_ref *= (nu0/sync_p_nu0)^(-3)
u_sync_ref *= (nu0/sync_p_nu0)^(-3)

i_dust_ref *= (nu0/dust_i_nu0)^(1.7)
q_dust_ref *= (nu0/dust_p_nu0)^(1.7)
u_dust_ref *= (nu0/dust_p_nu0)^(1.7)

xmap = dblarr(npix,3)
xmap[*,0] = i_dust_ref
xmap[*,1] = q_dust_ref
xmap[*,2] = u_dust_ref
ispice, xmap, cover, l_out, clt_dust, cle_dust, clb_dust, clte_dust

xmap = dblarr(npix,3)
xmap[*,0] = i_sync_ref
xmap[*,1] = q_sync_ref
xmap[*,2] = u_sync_ref
ispice, xmap, cover, l_out, clt_sync, cle_sync, clb_sync, clte_sync

;; wl = where( l_out ge 20 and l_out le nside)
wl = where( l_out ge 2 and l_out le 20)
dust_fit_t  = linfit( alog(l_out[wl]), alog( abs(clt_dust[wl])))
dust_fit_e  = linfit( alog(l_out[wl]), alog( abs(cle_dust[wl])))
dust_fit_b  = linfit( alog(l_out[wl]), alog( abs(clb_dust[wl])))
dust_fit_te = linfit( alog(l_out[wl]), alog( abs(clte_dust[wl])))

sync_fit_t  = linfit( alog(l_out[wl]), alog( abs(clt_sync[wl])))
sync_fit_e  = linfit( alog(l_out[wl]), alog( abs(cle_sync[wl])))
sync_fit_b  = linfit( alog(l_out[wl]), alog( abs(clb_sync[wl])))
sync_fit_te = linfit( alog(l_out[wl]), alog( abs(clte_sync[wl])))

col_e = 70
col_te = 150
col_b = 250
wind, 1, 1, /free, /large
my_multiplot, 2, 1, pp, pp1, /rev
plot_oo, l_out, l_out*(l_out+1)/(2*!dpi)*abs(clt_dust), xra=[1, 1000], yra=[1d-4, 1d4], position=pp1[0,*]
oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(cle_dust), col=col_e
oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clb_dust), col=col_b
oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clte_dust), col=col_te
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_t[0])*l_out^dust_fit_t[1]
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_e[0])*l_out^dust_fit_e[1], col=col_e
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_b[0])*l_out^dust_fit_b[1], col=col_b
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_te[0])*l_out^dust_fit_te[1], col=col_te
legendastro, 'Dust'

plot_oo, l_out, l_out*(l_out+1)/(2*!dpi)*abs(clt_sync), xra=[1, 1000], yra=[1d-10, 1d2], position=pp1[1,*], /noerase
oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(cle_sync), col=col_e
oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clb_sync), col=col_b
oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clte_sync), col=col_te
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_t[0])*l_out^sync_fit_t[1]
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_e[0])*l_out^sync_fit_e[1], col=col_e
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_b[0])*l_out^sync_fit_b[1], col=col_b
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_te[0])*l_out^sync_fit_te[1], col=col_te
legendastro, 'Sync'
stop

;; dust at nu0
clt  = exp(dust_fit_t[ 0])*l^dust_fit_t[ 1]
cle  = exp(dust_fit_e[ 0])*l^dust_fit_e[ 1]
clb  = exp(dust_fit_b[ 0])*l^dust_fit_b[ 1]
clte = exp(dust_fit_te[0])*l^dust_fit_te[1]
cls2mapsiqu, l, clt, cle, clb, clte, nx, res_arcmin/60., $
             dust_t, dust_q, dust_u, cu_t, cu_e, cu_b, cu_te

;; sync at nu0
clt  = exp(sync_fit_t[ 0])*l^sync_fit_t[ 1]
cle  = exp(sync_fit_e[ 0])*l^sync_fit_e[ 1]
clb  = exp(sync_fit_b[ 0])*l^sync_fit_b[ 1]
clte = exp(sync_fit_te[0])*l^sync_fit_te[1]

clte = clte < cle
clt  = clt > clte

cls2mapsiqu, l, clt, cle, clb, clte, nx, res_arcmin/60., $
             sync_t, sync_q, sync_u, cu_t, cu_e, cu_b, cu_te

wind, 1, 1, /free, /large
dp = {noerase:1}
my_multiplot, 3, 3, pp, pp1, /rev
imview, cmb_t, dp=dp, position=pp1[0,*], title='CMB T'
imview, cmb_q, dp=dp, position=pp1[1,*], title='CMB Q'
imview, cmb_u, dp=dp, position=pp1[2,*], title='CMB U'

imview, dust_t, dp=dp, position=pp1[3,*], title='DUST T'
imview, dust_q, dp=dp, position=pp1[4,*], title='DUST Q'
imview, dust_u, dp=dp, position=pp1[5,*], title='DUST U'

imview, sync_t, dp=dp, position=pp1[6,*], title='SYNC T'
imview, sync_q, dp=dp, position=pp1[7,*], title='SYNC Q'
imview, sync_u, dp=dp, position=pp1[8,*], title='SYNC U'
stop

;; ;;------------------------------------------------------------------------
;; ;; IF ONLY ONE COMPONENT
;; ;;-----------------------------------------------------------------------
;; ;; add non linearity, CMB only and acount for multiple detector observations
;; ata = dblarr(3,3)
;; atd = dblarr( long(nx)*long(ny), 3)
;; for idet=0, ndet-1 do begin
;;    cos2alpha = cos(2*alpha[idet])
;;    sin2alpha = sin(2*alpha[idet])
;; 
;;    ata[0,0] += 1.d0
;;    ata[1,0] += cos2alpha
;;    ata[2,0] += sin2alpha
;;    ata[1,1] += cos2alpha^2
;;    ata[2,1] += cos2alpha*sin2alpha
;;    ata[2,2] += sin2alpha^2
;; 
;;    ;; measure
;;    m = map_t + cos2alpha*map_q + sin2alpha*map_u
;; 
;;    m += epsilon*m^2
;; 
;;    atd[*,0] += m
;;    atd[*,1] += cos2alpha*m
;;    atd[*,2] += sin2alpha*m
;; endfor
;; ata[0,1] = ata[1,0]
;; ata[0,2] = ata[2,0]
;; ata[1,2] = ata[2,1]
;; 
;; atam1 = invert(ata)
;; 
;; ;; check
;; map_t = reform( (atam1##atd)[*,0], nx, ny)
;; map_q = reform( (atam1##atd)[*,1], nx, ny)
;; map_u = reform( (atam1##atd)[*,2], nx, ny)

;;-----------------------------------------------------------------------------------------------
;; Compute power spectra
s = size(map_q)
n        = s[1]
lmap_rad = float(n*res_arcmin*!arcmin2rad)

amn_q = fft( map_q, /double)
amn_u = fft( map_u, /double)
amn_e = amn_q*0.d0
amn_b = amn_q*0.d0

ic = complex( 0.00d0, 1.00d0)

for im=0L, n-1 do begin
   if im le n/2 then m1 = float(im) else m1 = float(n-im)
   for in=0L, n-1 do begin
      if in le n/2 then n1 = float(in) else n1 = float(n-in)
      um = m1/float(Lmap_rad)
      un = n1/float(Lmap_rad)
      u2 = um^2 + un^2
      if u2 eq 0 then begin
         amn_e[im,in] = 0.d0
         amn_b[im,in] = 0.d0
      endif else begin
         ;amn_q[im,in] = amn_e[im,in]*(um^2 - un^2)/u2 - amn_b[im,in]*2.0d0*um*un/u2
         ;amn_u[im,in] = amn_e[im,in]*2.0d0*um*un/u2   + amn_b[im,in]*(um^2 - un^2)/u2

         amn_e[im,in] =  amn_q[im,in]*(um^2 - un^2)/u2 + amn_u[im,in]*2.0d0*um*un/u2
         amn_b[im,in] = -amn_q[im,in]*2.0d0*um*un/u2   + amn_u[im,in]*(um^2 - un^2)/u2
      endelse
   endfor
endfor

map_e = double( fft( amn_e, /inverse, /double))
map_b = double( fft( amn_b, /inverse, /double))

delta_l_over_l = 0.1 ; default
delta_l_over_l = 0.05
ipoker, map_t, res_arcmin, k, pk_i, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, map_e, res_arcmin, k, pk_e, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, map_b, res_arcmin, k, pk_B, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, map_t, res_arcmin, map1=map_e, k, pk_te, /bypass, /rem, delta_l_over_l=delta_l_over_l

wind, 1, 1, /free, /large
my_multiplot, 3, 2, pp, pp1, /rev, gap_x=0.05
chars = 0.7
dp = {noerase:1, charsize:chars, charbar:chars}
imview, map_t, dp=dp, position=pp1[0,*], title='T'
imview, map_q, dp=dp, position=pp1[1,*], title='Q'
imview, map_u, dp=dp, position=pp1[2,*], title='U'

imview, map_e, dp=dp, position=pp1[4,*], title='E'
imview, map_b, dp=dp, position=pp1[5,*], title='B'

yra = [1d-6,1d4] ; cmb
yra = [1d-10, 1d3] ; place holder foregrounds
xra = [10,2000]
psym = 8 ; -8
syms = 0.5
col_e = 70
col_b = 250
col_te = 150
stop
plot_oo, k, k*(k+1)/(2*!dpi)*pk_i, yra=yra, /ys, xra=xra, /xs, $
         psym=psym, syms=syms, $
         position=pp[0,1,*], /noerase, $
         xtitle='Multipole l', ytitle='l(l+1)/2!7p!3 C!dl'
oplot, l, l*(l+1)/(2*!dpi)*clt
oplot, k, k*(k+1)/(2*!dpi)*abs(pk_te), col=col_te, psym=psym, syms=syms
oplot, l, l*(l+1)/(2*!dpi)*abs(clte), col=col_te
oplot, l, l*(l+1)/(2*!dpi)*clb_in, col=0
oplot, k, k*(k+1)/(2*!dpi)*pk_e, col=col_e, psym=psym, syms=syms
oplot, l, l*(l+1)/(2*!dpi)*cle, col=col_e
oplot, k, k*(k+1)/(2*!dpi)*pk_b, col=col_b, psym=psym, syms=syms
oplot, l, l*(l+1)/(2*!dpi)*clb, col=col_b
legendastro, 'Epsilon = '+strtrim(epsilon,2)
legendastro, ['TT', 'EE', 'TE', 'BB'], col=[!p.color, col_e, col_te, col_b], /bottom, line=0
;; 
;; 
;; 
;; for imc=0, nmc-1 do begin
;;    isynfast, input_cl_file, cmb_maps, nlmax=3*nside-1, nside=nside
;;    i_cmb = reform( cmb_maps[*,0])
;;    q_cmb = reform( cmb_maps[*,1])
;;    u_cmb = reform( cmb_maps[*,2])
;; 
;;    ata = dblarr(9,9)
;;    atd = dblarr(npix,9)
;;    
;;    for inu=0, n_nu-1 do begin
;;       rj2k = rj2thermo(nu[inu])
;;       
;;       for ialpha=0, ndet-1 do begin
;;          cos2alpha = cos(2.*alpha[ialpha])
;;          sin2alpha = sin(2.*alpha[ialpha])
;;          
;;          ;; Fill ata
;;          ata[0,0] += 1.d0
;;          ata[1,0] += cos2alpha
;;          ata[2,0] += sin2alpha
;;          ata[3,0] += rj2k * (nu[inu]/dust_i_nu0)^beta_dust
;;          ata[4,0] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha
;;          ata[5,0] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha
;;          ata[6,0] += rj2k * (nu[inu]/sync_i_nu0)^beta_sync
;;          ata[7,0] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha
;;          ata[8,0] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha
;; 
;;          ata[1,1] += cos2alpha^2
;;          ata[2,1] += cos2alpha*sin2alpha
;;          ata[3,1] += rj2k * (nu[inu]/dust_i_nu0)^beta_dust*cos2alpha
;;          ata[4,1] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha^2
;;          ata[5,1] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha*cos2alpha
;;          ata[6,1] += rj2k * (nu[inu]/sync_i_nu0)^beta_sync*cos2alpha
;;          ata[7,1] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha^2
;;          ata[8,1] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha*cos2alpha
;; 
;;          ata[2,2] += sin2alpha^2
;;          ata[3,2] += rj2k * (nu[inu]/dust_i_nu0)^beta_dust*sin2alpha
;;          ata[4,2] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha*sin2alpha
;;          ata[5,2] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha^2
;;          ata[6,2] += rj2k * (nu[inu]/sync_i_nu0)^beta_sync*sin2alpha
;;          ata[7,2] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha*sin2alpha
;;          ata[8,2] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha^2
;; 
;;          ata[3,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^(2.*beta_dust)
;;          ata[4,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha
;;          ata[5,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha
;;          ata[6,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/sync_i_nu0)^beta_sync
;;          ata[7,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha
;;          ata[8,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha
;; 
;;          ata[4,4] += rj2k^2 * (nu[inu]/dust_p_nu0)^(2.*beta_dust)*cos2alpha^2
;;          ata[5,4] += rj2k^2 * (nu[inu]/dust_p_nu0)^(2.*beta_dust)*cos2alpha*sin2alpha
;;          ata[6,4] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_i_nu0)^beta_sync*cos2alpha
;;          ata[7,4] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha^2
;;          ata[8,4] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha*sin2alpha
;; 
;;          ata[5,5] += rj2k^2 * (nu[inu]/dust_p_nu0)^(2.*beta_dust)*sin2alpha^2
;;          ata[6,5] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_i_nu0)^beta_sync*sin2alpha
;;          ata[7,5] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha*sin2alpha
;;          ata[8,5] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha^2
;; 
;;          ata[6,6] += rj2k^2 * (nu[inu]/sync_i_nu0)^(2.*beta_sync)
;;          ata[7,6] += rj2k^2 * (nu[inu]/sync_i_nu0)^beta_sync * (nu[inu]/sync_p_nu0)^beta_sync* cos2alpha
;;          ata[8,6] += rj2k^2 * (nu[inu]/sync_i_nu0)^beta_sync * (nu[inu]/sync_p_nu0)^beta_sync* sin2alpha
;; 
;;          ata[7,7] += rj2k^2 * (nu[inu]/sync_p_nu0)^(2.*beta_sync)*cos2alpha^2
;;          ata[8,7] += rj2k^2 * (nu[inu]/sync_p_nu0)^(2.*beta_sync)*cos2alpha*sin2alpha
;; 
;;          ata[8,8] += rj2k^2 * (nu[inu]/sync_p_nu0)^(2.*beta_sync)*sin2alpha^2
;; 
;;          ;;------------------------------
;;          ;; Define A^tD
;;          
;; ;;         ;; Scale components to the current frequency in microK_RJ then to microK_CMB
;; ;;         i_dust = rj2k*(nu[inu]/dust_i_nu0)^beta_dust*i_dust_ref
;; ;;         q_dust = rj2k*(nu[inu]/dust_p_nu0)^beta_dust*q_dust_ref
;; ;;         u_dust = rj2k*(nu[inu]/dust_p_nu0)^beta_dust*u_dust_ref
;; ;;         i_sync = rj2k*(nu[inu]/sync_i_nu0)^beta_sync*i_sync_ref
;; ;;         q_sync = rj2k*(nu[inu]/sync_p_nu0)^beta_sync*q_sync_ref
;; ;;         u_sync = rj2k*(nu[inu]/sync_p_nu0)^beta_sync*u_sync_ref
;; 
;;          ;; Scale components to the current frequency in microK_RJ
;;          i_dust = (nu[inu]/dust_i_nu0)^beta_dust*i_dust_ref
;;          q_dust = (nu[inu]/dust_p_nu0)^beta_dust*q_dust_ref
;;          u_dust = (nu[inu]/dust_p_nu0)^beta_dust*u_dust_ref
;;          i_sync = (nu[inu]/sync_i_nu0)^beta_sync*i_sync_ref
;;          q_sync = (nu[inu]/sync_p_nu0)^beta_sync*q_sync_ref
;;          u_sync = (nu[inu]/sync_p_nu0)^beta_sync*u_sync_ref
;; 
;;          ;; Convert to Mjy/sr
;;          lambda_microns = !const.c/(nu[inu]*1d9)*1d6
;;          convert_millik_megajy, lambda_microns, i_dust/1000., i_dust_megajy, /rj
;;          convert_millik_megajy, lambda_microns, q_dust/1000., q_dust_megajy, /rj
;;          convert_millik_megajy, lambda_microns, u_dust/1000., u_dust_megajy, /rj
;;          convert_millik_megajy, lambda_microns, i_sync/1000., i_sync_megajy, /rj
;;          convert_millik_megajy, lambda_microns, q_sync/1000., q_sync_megajy, /rj
;;          convert_millik_megajy, lambda_microns, u_sync/1000., u_sync_megajy, /rj
;; 
;;          convert_millik_megajy, lambda_microns, i_cmb/1000., i_cmb_megajy, /cmb
;;          convert_millik_megajy, lambda_microns, q_cmb/1000., q_cmb_megajy, /cmb
;;          convert_millik_megajy, lambda_microns, u_cmb/1000., u_cmb_megajy, /cmb
;; 
;;          ;; Derive measured flux in Mjy
;;          fi = (i_cmb_megajy + i_dust_megajy + i_sync_megajy)
;;          fq = (q_cmb_megajy + q_dust_megajy + q_sync_megajy)
;;          fu = (u_cmb_megajy + u_dust_megajy + u_sync_megajy)
;; 
;;          ;; Account for non linearity
;;          m = fi + cos2alpha*fq + sin2alpha*fu +$
;;              epsilon * (fi + cos2alpha*fq + sin2alpha*fu)^2
;; 
;;          ;; Convert back to thermodynamic temperature
;;          convert_megajy_millik, lambda_microns, m, m_mk, /cmb
;;          m_microK = m_mk*1d3
;; 
;;          atd[*,0] += m_microK
;;          atd[*,1] += m_microK*cos2alpha
;;          atd[*,2] += m_microK*sin2alpha
;;          atd[*,3] += rj2k * m_microK*(nu[inu]/dust_i_nu0)^beta_dust
;;          atd[*,4] += rj2k * m_microK*(nu[inu]/dust_p_nu0)^beta_dust*cos2alpha
;;          atd[*,5] += rj2k * m_microK*(nu[inu]/dust_p_nu0)^beta_dust*sin2alpha
;;          atd[*,6] += rj2k * m_microK*(nu[inu]/sync_i_nu0)^beta_sync
;;          atd[*,7] += rj2k * m_microK*(nu[inu]/sync_p_nu0)^beta_sync*cos2alpha
;;          atd[*,8] += rj2k * m_microK*(nu[inu]/sync_p_nu0)^beta_sync*sin2alpha
;;          
;;       endfor
;;    endfor
;; ;; complete ata
;;    for i=1, 8 do begin
;;       for j=0, i-1 do begin
;;          ata[j,i] = ata[i,j]
;;       endfor
;;    endfor
;; 
;;    s = invert(ata)##atd
;;    out_cmb_i  = reform( s[*,0])
;;    out_cmb_q  = reform( s[*,1])
;;    out_cmb_u  = reform( s[*,2])
;; ;;    out_dust_i = reform( s[*,3])
;; ;;    out_dust_q = reform( s[*,4])
;; ;;    out_dust_u = reform( s[*,5])
;; ;;    out_sync_i = reform( s[*,6])
;; ;;    out_sync_q = reform( s[*,7])
;; ;;    out_sync_u = reform( s[*,8])
;; ;; 
;; 
;; ;; mollview, cover*(out_cmb_i-i_cmb), title='CMB I out-in, epsilon='+strtrim(epsilon,2)
;; ;; mollview, cover*(out_cmb_q-q_cmb), title='CMB Q out-in, epsilon='+strtrim(epsilon,2)
;; ;;   stop
;; ;; mollview, out_cmb_q-q_cmb, title='CMB Q out-in'
;; ;; mollview, cover*(out_dust_i-i_dust_ref), title='Dust (ref nu) out-in'
;; ;; mollview, out_sync_q-q_sync_ref, title='Sync (ref nu) out-in (RJ)'
;; ;; print, "HERE"
;; ;; stop
;; 
;;    ;; Power spectra
;;    xmap = dblarr(npix,3)
;;    xmap[*,0] = out_cmb_i
;;    xmap[*,1] = out_cmb_q
;;    xmap[*,2] = out_cmb_u
;;    ispice, xmap, cover, l_out, clt_out, cle_out, clb_out, clte_out
;;    l_out    = l_out[   2:*]
;;    clt_out  = clt_out[ 2:*]
;;    cle_out  = cle_out[ 2:*]
;;    clb_out  = clb_out[ 2:*]
;;    clte_out = clte_out[2:*]
;;    
;; ;; Correct for Healpix's pixel window function
;;    wl = mrdfits( !healpix.path.data+"/pixel_window_n"+string(nside,form='(I4.4)')+".fits",1)
;;    clt_out  = clt_out[ 2:*]/wl.temperature^2
;;    cle_out  = cle_out[ 2:*]/wl.polarization^2
;;    clb_out  = clb_out[ 2:*]/wl.polarization^2
;;    clte_out = clte_out[2:*]/wl.temperature/wl.polarization
;; 
;;    if imc eq 0 then begin
;;       clt_out_res  = dblarr(nmc,n_elements(clt_out))
;;       cle_out_res  = clt_out_res
;;       clb_out_res  = clt_out_res
;;       clte_out_res = clt_out_res
;;    endif
;;    clt_out_res[imc,*]  = clt_out
;;    cle_out_res[imc,*]  = cle_out
;;    clb_out_res[imc,*]  = clb_out
;;    clte_out_res[imc,*] = clte_out
;; 
;; ;;    ;; One realization
;; ;;    fl = l_out*(l_out+1)/(2.*!dpi)
;; ;;    xra = [1,3*nside]
;; ;;    yra = [1d-8,1d4]
;; ;;    col_ee = 70
;; ;;    col_bb = 250
;; ;;    col_te = 150
;; ;;    wind, 1, 1, /free, /large
;; ;;    plot, xra, yra, xtitle='!12l!3', ytitle='!12l(l+1)C!dl!n/2!7p!3', $
;; ;;          /nodata, /xlog, /ylog, /xs, /ys
;; ;;    legendastro, 'Epsilon '+strtrim(epsilon,2), /bottom
;; ;;    oplot, l, clt
;; ;;    oplot, l, cle
;; ;;    oplot, l, clb
;; ;;    oplot, l, abs(clte)
;; ;;    oplot, l_out, fl*clt_out
;; ;;    oplot, l_out, fl*cle_out, col=col_ee
;; ;;    oplot, l_out, fl*clb_out, col=col_bb
;; ;;    oplot, l_out, fl*abs(clte_out), col=col_te
;; ;;    legendastro, ['T', 'E', 'B', 'TE'], col=[!p.color, col_ee, col_bb, col_te], line=0
;; ;;    stop
;; endfor
;; 
;; ;; Plot average results
;; mc_reduce, clt_out_res, clt_out_avg, sigma_clt_out_avg
;; mc_reduce, cle_out_res, cle_out_avg, sigma_cle_out_avg
;; mc_reduce, clb_out_res, clb_out_avg, sigma_clb_out_avg
;; mc_reduce, clte_out_res, clte_out_avg, sigma_clte_out_avg
;; 
;; ;; Derive constraint on epsilon
;; l_ref = 100
;; w_out = (where( abs(l_out-l_ref) eq min( abs(l_out-l_ref))))[0]
;; w     = (where( abs(l-l_ref) eq min( abs(l-l_ref))))[0]
;; print, "epsilon, clb_out_avg(l="+strtrim(l_ref,2)+")/cl_cmb = ", $
;;        clb_out_avg[w]/(clb[w]*2*!dpi/(l[w]*(l[w]+1)))
;; 
;; ;; Final plot
;; fl = l_out*(l_out+1)/(2.*!dpi)
;; xra = [1,3*nside]
;; yra = [1d-8,1d4]
;; col_ee = 70
;; col_bb = 250
;; col_te = 150
;; wind, 1, 1, /free, /large
;; outplot, file='non_linear_fg_residuals_cl_epsilon_'+strtrim(epsilon,2), png=png, ps=ps
;; plot, xra, yra, xtitle='!12l!3', ytitle='!12l(l+1)C!dl!n/2!7p!3', $
;;       /nodata, /xlog, /ylog, /xs, /ys
;; legendastro, 'Epsilon '+strtrim(epsilon,2), /bottom
;; oplot, l, clt
;; oplot, l, cle
;; oplot, l, clb
;; oplot, l, abs(clte)
;; oplot, l_out, fl*clt_out_avg
;; oplot, l_out, fl*cle_out_avg, col=col_ee
;; oplot, l_out, fl*clb_out_avg, col=col_bb
;; oplot, l_out, fl*abs(clte_out_avg), col=col_te
;; legendastro, ['T', 'E', 'B', 'TE'], col=[!p.color, col_ee, col_bb, col_te], line=0
;; legendastro, 'Nmc = '+strtrim(nmc,2), /right
;; outplot, /close, /verb


end
