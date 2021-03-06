
;; Simple component separation to assess impact of non linearity

;; Assumed non-linearity parameter
epsilon = 5d-2 ; 0.d0 ; 1d-4 ;0.d0 ;1d-3 ;0.5d-1 ;1d-6 ;0.d0
epsilon = [5d-2, 1d-1] 
epsilon = [5d-3, 5d-2, 1d-1] ;plot these

eps_min = 5d-3
eps_max = 0.20  ;0.25
neps    = 15
epsilon = dindgen(neps)*(eps_max-eps_min)/(neps-1) + eps_min

;;epsilon = [1d-3, 3d-3, 5d-3, 1d-2, 3d-2, 5d-2, 1d-1, 3d-1, 5d-1, 1, 1.3, 1.5, 2d0] 
neps    = n_elements(epsilon)

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
r = 8 ; 9
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

;; Number of components
ncomp = 3

;; Number of polarization states
nstokes = 3


;; MC realizations
nmc = 5 ; 30 ;10 ; 50 ;30

;;============================================
n_nu  = n_elements(nu)
alpha = alpha_det_deg*!dtor
ndet  = n_elements(alpha)

;; Quick check
readcol, '$SK_DIR/Cl/cmb_totCls_r0.001.dat', l, clt, cle, clb, clte, $
         format='D,D,D,D', comment='#'
;; prefactor included in the CAMB's output files
fl = l*(l+1)/(2*!dpi)
cmb_clt  = clt/fl
cmb_cle  = cle/fl
cmb_clb_in = clb/fl
cmb_clte = clte/fl
;; No B in the simulations
cmb_clb = cle*0.d0
cls2mapsiqu, l, cmb_clt, cmb_cle, cmb_clb, cmb_clte, nx, res_arcmin/60., $
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
;;rj2k_nu0 = rj2thermo(nu0)

beta_sync = -3.
i_sync = i_sync_ref * (nu0/sync_i_nu0)^beta_sync;; * rj2k_nu0
q_sync = q_sync_ref * (nu0/sync_p_nu0)^beta_sync;; * rj2k_nu0
u_sync = u_sync_ref * (nu0/sync_p_nu0)^beta_sync;; * rj2k_nu0

beta_dust = 1.5
i_dust = i_dust_ref * (nu0/dust_i_nu0)^beta_dust;;  * rj2k_nu0
q_dust = q_dust_ref * (nu0/dust_p_nu0)^beta_dust;;  * rj2k_nu0
u_dust = u_dust_ref * (nu0/dust_p_nu0)^beta_dust;;  * rj2k_nu0

xmap = dblarr(npix,3)
xmap[*,0] = i_dust
xmap[*,1] = q_dust
xmap[*,2] = u_dust
ispice, xmap, cover, l_out, clt_dust, cle_dust, clb_dust, clte_dust

xmap = dblarr(npix,3)
xmap[*,0] = i_sync
xmap[*,1] = q_sync
xmap[*,2] = u_sync
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
;; wind, 1, 1, /free, /large
;; my_multiplot, 2, 1, pp, pp1, /rev
;; plot_oo, l_out, l_out*(l_out+1)/(2*!dpi)*abs(clt_dust), xra=[1, 1000], yra=[1d-4, 1d4], position=pp1[0,*]
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(cle_dust), col=col_e
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clb_dust), col=col_b
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clte_dust), col=col_te
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_t[0])*l_out^dust_fit_t[1]
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_e[0])*l_out^dust_fit_e[1], col=col_e
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_b[0])*l_out^dust_fit_b[1], col=col_b
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_te[0])*l_out^dust_fit_te[1], col=col_te
;; legendastro, 'Dust'

;; plot_oo, l_out, l_out*(l_out+1)/(2*!dpi)*abs(clt_sync), xra=[1, 1000], yra=[1d-10, 1d2], position=pp1[1,*], /noerase
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(cle_sync), col=col_e
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clb_sync), col=col_b
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clte_sync), col=col_te
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_t[0])*l_out^sync_fit_t[1]
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_e[0])*l_out^sync_fit_e[1], col=col_e
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_b[0])*l_out^sync_fit_b[1], col=col_b
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_te[0])*l_out^sync_fit_te[1], col=col_te
;; legendastro, 'Sync'
;stop

;; dust at nu0
dust_clt  = exp(dust_fit_t[ 0])*l^dust_fit_t[ 1]
dust_cle  = exp(dust_fit_e[ 0])*l^dust_fit_e[ 1]
dust_clb  = exp(dust_fit_b[ 0])*l^dust_fit_b[ 1]
dust_clte = exp(dust_fit_te[0])*l^dust_fit_te[1]
cls2mapsiqu, l, dust_clt, dust_cle, dust_clb, dust_clte, nx, res_arcmin/60., $
             dust_t, dust_q, dust_u, cu_t, cu_e, cu_b, cu_te

;; sync at nu0
sync_clt  = exp(sync_fit_t[ 0])*l^sync_fit_t[ 1]
sync_cle  = exp(sync_fit_e[ 0])*l^sync_fit_e[ 1]
sync_clb  = exp(sync_fit_b[ 0])*l^sync_fit_b[ 1]
sync_clte = exp(sync_fit_te[0])*l^sync_fit_te[1]
sync_clte = sync_clte < sync_cle
sync_clt  = sync_clt  > sync_clte
cls2mapsiqu, l, sync_clt, sync_cle, sync_clb, sync_clte, nx, res_arcmin/60., $
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

;; Xcheck angular power spectra of simulated maps
s        = size(cmb_t)
n        = s[1]
lmap_rad = float(n*res_arcmin*!arcmin2rad)

qu2eb, cmb_q, cmb_u,   res_arcmin, cmb_e, cmb_b
qu2eb, dust_q, dust_u, res_arcmin, dust_e, dust_b
qu2eb, sync_q, sync_u, res_arcmin, sync_e, sync_b

delta_l_over_l = 0.05
ipoker, cmb_t, res_arcmin, k, pk_cmb_t, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, cmb_e, res_arcmin, k, pk_cmb_e, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, cmb_b, res_arcmin, k, pk_cmb_b, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, cmb_t, res_arcmin, map1=cmb_e, k, pk_cmb_te, /bypass, /rem, delta_l_over_l=delta_l_over_l

ipoker, dust_t, res_arcmin, k, pk_dust_t, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, dust_e, res_arcmin, k, pk_dust_e, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, dust_b, res_arcmin, k, pk_dust_b, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, dust_t, res_arcmin, map1=dust_e, k, pk_dust_te, /bypass, /rem, delta_l_over_l=delta_l_over_l

ipoker, sync_t, res_arcmin, k, pk_sync_t, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, sync_e, res_arcmin, k, pk_sync_e, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, sync_b, res_arcmin, k, pk_sync_b, /bypass, /rem, delta_l_over_l=delta_l_over_l
ipoker, sync_t, res_arcmin, map1=sync_e, k, pk_sync_te, /bypass, /rem, delta_l_over_l=delta_l_over_l

;; wind, 1, 1, /free, /large
;; my_multiplot, 1, 3, pp, pp1, /rev
;; chars = 0.7
;; dp = {noerase:1, charsize:chars, charbar:chars}
;; yra = array2range(k*(k+1)/(2*!dpi)*[pk_cmb_t,pk_cmb_b])
;; psym = 8 ; -8
;; syms = 0.5
;; col_e = 70
;; col_b = 250
;; col_te = 150
;; plot_oo, k, k*(k+1)/(2*!dpi)*pk_cmb_t, yra=yra, /ys, xra=xra, /xs, $
;;          psym=psym, syms=syms, $
;;          position=pp1[0,*], /noerase, $
;;          xtitle='Multipole l', ytitle='l(l+1)/2!7p!3 C!dl'
;; oplot, l, l*(l+1)/(2*!dpi)*cmb_clt
;; oplot, k, k*(k+1)/(2*!dpi)*abs(pk_cmb_te), col=col_te, psym=psym, syms=syms
;; oplot, l, l*(l+1)/(2*!dpi)*abs(cmb_clte), col=col_te
;; oplot, l, l*(l+1)/(2*!dpi)*cmb_clb, col=0
;; oplot, k, k*(k+1)/(2*!dpi)*pk_cmb_e, col=col_e, psym=psym, syms=syms
;; oplot, l, l*(l+1)/(2*!dpi)*cmb_cle, col=col_e
;; oplot, k, k*(k+1)/(2*!dpi)*pk_cmb_b, col=col_b, psym=psym, syms=syms
;; oplot, l, l*(l+1)/(2*!dpi)*cmb_clb, col=col_b
;; legendastro, 'Epsilon = '+strtrim(epsilon,2)
;; legendastro, ['TT', 'EE', 'TE', 'BB'], col=[!p.color, col_e, col_te, col_b], /bottom, line=0

;; yra = array2range(k*(k+1)/(2*!dpi)*[pk_dust_t,pk_dust_b])
;; plot_oo, k, k*(k+1)/(2*!dpi)*pk_dust_t, yra=yra, /ys, xra=xra, /xs, $
;;          psym=psym, syms=syms, $
;;          position=pp1[1,*], /noerase, $
;;          xtitle='Multipole l', ytitle='l(l+1)/2!7p!3 C!dl'
;; oplot, l, l*(l+1)/(2*!dpi)*dust_clt
;; oplot, k, k*(k+1)/(2*!dpi)*abs(pk_dust_te), col=col_te, psym=psym, syms=syms
;; oplot, l, l*(l+1)/(2*!dpi)*abs(dust_clte), col=col_te
;; oplot, l, l*(l+1)/(2*!dpi)*dust_clb, col=0
;; oplot, k, k*(k+1)/(2*!dpi)*pk_dust_e, col=col_e, psym=psym, syms=syms
;; oplot, l, l*(l+1)/(2*!dpi)*dust_cle, col=col_e
;; oplot, k, k*(k+1)/(2*!dpi)*pk_dust_b, col=col_b, psym=psym, syms=syms
;; oplot, l, l*(l+1)/(2*!dpi)*dust_clb, col=col_b
;; legendastro, 'Epsilon = '+strtrim(epsilon,2)
;; legendastro, ['TT', 'EE', 'TE', 'BB'], col=[!p.color, col_e, col_te, col_b], /bottom, line=0

;; yra = array2range(k*(k+1)/(2*!dpi)*[pk_sync_t,pk_sync_b])
;; plot_oo, k, k*(k+1)/(2*!dpi)*pk_sync_t, yra=yra, /ys, xra=xra, /xs, $
;;          psym=psym, syms=syms, $
;;          position=pp1[2,*], /noerase, $
;;          xtitle='Multipole l', ytitle='l(l+1)/2!7p!3 C!dl'
;; oplot, l, l*(l+1)/(2*!dpi)*sync_clt
;; oplot, k, k*(k+1)/(2*!dpi)*abs(pk_sync_te), col=col_te, psym=psym, syms=syms
;; oplot, l, l*(l+1)/(2*!dpi)*abs(sync_clte), col=col_te
;; oplot, l, l*(l+1)/(2*!dpi)*sync_clb, col=0
;; oplot, k, k*(k+1)/(2*!dpi)*pk_sync_e, col=col_e, psym=psym, syms=syms
;; oplot, l, l*(l+1)/(2*!dpi)*sync_cle, col=col_e
;; oplot, k, k*(k+1)/(2*!dpi)*pk_sync_b, col=col_b, psym=psym, syms=syms
;; oplot, l, l*(l+1)/(2*!dpi)*sync_clb, col=col_b
;; legendastro, 'Epsilon = '+strtrim(epsilon,2)
;; legendastro, ['TT', 'EE', 'TE', 'BB'], col=[!p.color, col_e, col_te, col_b], /bottom, line=0
;stop


;;=================================== Main loop =======================

ata = dblarr(ncomp*nstokes, ncomp*nstokes)
;atd = dblarr(ncomp*nstokes)
out_cmb_t = dblarr(nx, ny, neps)
out_cmb_q = dblarr(nx, ny, neps)
out_cmb_u = dblarr(nx, ny, neps)

for inu=0, n_elements(nu)-1 do begin
   rj2k = rj2thermo(nu[inu])
   for ialpha=0, ndet-1 do begin
      cos2alpha = cos(2*alpha[ialpha])
      sin2alpha = sin(2*alpha[ialpha])

      ata[0,0] += 1.d0
      ata[1,0] += cos2alpha
      ata[2,0] += sin2alpha
      ata[3,0] += rj2k * (nu[inu]/nu0)^beta_dust
      ata[4,0] += rj2k * (nu[inu]/nu0)^beta_dust*cos2alpha
      ata[5,0] += rj2k * (nu[inu]/nu0)^beta_dust*sin2alpha
      ata[6,0] += rj2k * (nu[inu]/nu0)^beta_sync
      ata[7,0] += rj2k * (nu[inu]/nu0)^beta_sync*cos2alpha
      ata[8,0] += rj2k * (nu[inu]/nu0)^beta_sync*sin2alpha

      ata[1,1] += cos2alpha^2
      ata[2,1] += cos2alpha*sin2alpha
      ata[3,1] += rj2k * (nu[inu]/nu0)^beta_dust*cos2alpha
      ata[4,1] += rj2k * (nu[inu]/nu0)^beta_dust*cos2alpha^2
      ata[5,1] += rj2k * (nu[inu]/nu0)^beta_dust*sin2alpha*cos2alpha
      ata[6,1] += rj2k * (nu[inu]/nu0)^beta_sync*cos2alpha
      ata[7,1] += rj2k * (nu[inu]/nu0)^beta_sync*cos2alpha^2
      ata[8,1] += rj2k * (nu[inu]/nu0)^beta_sync*sin2alpha*cos2alpha

      ata[2,2] += sin2alpha^2
      ata[3,2] += rj2k * (nu[inu]/nu0)^beta_dust*sin2alpha
      ata[4,2] += rj2k * (nu[inu]/nu0)^beta_dust*cos2alpha*sin2alpha
      ata[5,2] += rj2k * (nu[inu]/nu0)^beta_dust*sin2alpha^2
      ata[6,2] += rj2k * (nu[inu]/nu0)^beta_sync*sin2alpha
      ata[7,2] += rj2k * (nu[inu]/nu0)^beta_sync*cos2alpha*sin2alpha
      ata[8,2] += rj2k * (nu[inu]/nu0)^beta_sync*sin2alpha^2

      ata[3,3] += rj2k^2 * (nu[inu]/nu0)^(2.*beta_dust)
      ata[4,3] += rj2k^2 * (nu[inu]/nu0)^beta_dust * (nu[inu]/nu0)^beta_dust*cos2alpha
      ata[5,3] += rj2k^2 * (nu[inu]/nu0)^beta_dust * (nu[inu]/nu0)^beta_dust*sin2alpha
      ata[6,3] += rj2k^2 * (nu[inu]/nu0)^beta_dust * (nu[inu]/nu0)^beta_sync
      ata[7,3] += rj2k^2 * (nu[inu]/nu0)^beta_dust * (nu[inu]/nu0)^beta_sync*cos2alpha
      ata[8,3] += rj2k^2 * (nu[inu]/nu0)^beta_dust * (nu[inu]/nu0)^beta_sync*sin2alpha

      ata[4,4] += rj2k^2 * (nu[inu]/nu0)^(2.*beta_dust)*cos2alpha^2
      ata[5,4] += rj2k^2 * (nu[inu]/nu0)^(2.*beta_dust)*cos2alpha*sin2alpha
      ata[6,4] += rj2k^2 * (nu[inu]/nu0)^beta_dust * (nu[inu]/nu0)^beta_sync*cos2alpha
      ata[7,4] += rj2k^2 * (nu[inu]/nu0)^beta_dust * (nu[inu]/nu0)^beta_sync*cos2alpha^2
      ata[8,4] += rj2k^2 * (nu[inu]/nu0)^beta_dust * (nu[inu]/nu0)^beta_sync*cos2alpha*sin2alpha
 
      ata[5,5] += rj2k^2 * (nu[inu]/nu0)^(2.*beta_dust)*sin2alpha^2
      ata[6,5] += rj2k^2 * (nu[inu]/nu0)^beta_dust * (nu[inu]/nu0)^beta_sync*sin2alpha
      ata[7,5] += rj2k^2 * (nu[inu]/nu0)^beta_dust * (nu[inu]/nu0)^beta_sync*cos2alpha*sin2alpha
      ata[8,5] += rj2k^2 * (nu[inu]/nu0)^beta_dust * (nu[inu]/nu0)^beta_sync*sin2alpha^2
      
      ata[6,6] += rj2k^2 * (nu[inu]/nu0)^(2.*beta_sync)
      ata[7,6] += rj2k^2 * (nu[inu]/nu0)^beta_sync * (nu[inu]/nu0)^beta_sync* cos2alpha
      ata[8,6] += rj2k^2 * (nu[inu]/nu0)^beta_sync * (nu[inu]/nu0)^beta_sync* sin2alpha
      
      ata[7,7] += rj2k^2 * (nu[inu]/nu0)^(2.*beta_sync)*cos2alpha^2
      ata[8,7] += rj2k^2 * (nu[inu]/nu0)^(2.*beta_sync)*cos2alpha*sin2alpha
      
      ata[8,8] += rj2k^2 * (nu[inu]/nu0)^(2.*beta_sync)*sin2alpha^2
   endfor
endfor
;; finish ata
for i=1, 8 do begin
   for j=0, i-1 do begin
      ata[j,i] = ata[i,j]
   endfor
endfor
;; invert
atam1 = invert(ata)      

nk = n_elements(k)
clt_out_res = dblarr(nmc,nk,neps)
cle_out_res = dblarr(nmc,nk,neps)
clb_out_res = dblarr(nmc,nk,neps)
clte_out_res = dblarr(nmc,nk,neps)
;m_eps        = dblarr(neps)
for imc=0, nmc-1 do begin
;   print, strtrim(imc,2)+"/"+strtrim(nmc-1,2)

   ;; cmb spectra are in microK_CMB^2
   cls2mapsiqu, l, cmb_clt, cmb_cle, cmb_clb, cmb_clte, nx, res_arcmin/60., $
                cmb_t, cmb_q, cmb_u

   ;; dust and sync spectra are in microK_RJ^2
   cls2mapsiqu, l, dust_clt, dust_cle, dust_clb, dust_clte, nx, res_arcmin/60., $
                nu0_dust_t, nu0_dust_q, nu0_dust_u
   cls2mapsiqu, l, sync_clt, sync_cle, sync_clb, sync_clte, nx, res_arcmin/60., $
                nu0_sync_t, nu0_sync_q, nu0_sync_u, cu_t, cu_e, cu_b, cu_te

   npix = n_elements(cmb_t)
   atd = dblarr(npix,ncomp*nstokes)
   atd = dblarr(npix,ncomp*nstokes,neps)
   s   = dblarr(npix,ncomp*nstokes,neps)

   for inu=0, n_elements(nu)-1 do begin
      rj2k = rj2thermo(nu[inu])
      
      dust_t = rj2k * (nu[inu]/nu0)^beta_dust * nu0_dust_t
      dust_q = rj2k * (nu[inu]/nu0)^beta_dust * nu0_dust_q
      dust_u = rj2k * (nu[inu]/nu0)^beta_dust * nu0_dust_u
      
      sync_t = rj2k * (nu[inu]/nu0)^beta_sync * nu0_sync_t
      sync_q = rj2k * (nu[inu]/nu0)^beta_sync * nu0_sync_q
      sync_u = rj2k * (nu[inu]/nu0)^beta_sync * nu0_sync_u
      
      for ialpha=0, ndet-1 do begin
         cos2alpha = cos(2*alpha[ialpha])
         sin2alpha = sin(2*alpha[ialpha])

         ;; Convert to Mjy/sr
         lambda_microns = !const.c/(nu[inu]*1d9)*1d6
         convert_millik_megajy, lambda_microns, cmb_t/1000., cmb_t_megajy, /cmb
         convert_millik_megajy, lambda_microns, cmb_q/1000., cmb_q_megajy, /cmb
         convert_millik_megajy, lambda_microns, cmb_u/1000., cmb_u_megajy, /cmb

         convert_millik_megajy, lambda_microns, dust_t/1000., dust_t_megajy, /cmb
         convert_millik_megajy, lambda_microns, dust_q/1000., dust_q_megajy, /cmb
         convert_millik_megajy, lambda_microns, dust_u/1000., dust_u_megajy, /cmb

         convert_millik_megajy, lambda_microns, sync_t/1000., sync_t_megajy, /cmb
         convert_millik_megajy, lambda_microns, sync_q/1000., sync_q_megajy, /cmb
         convert_millik_megajy, lambda_microns, sync_u/1000., sync_u_megajy, /cmb

         ;; Effective Stokes parameters
         flux_i = cmb_t_megajy + dust_t_megajy + sync_t_megajy
         flux_q = cmb_q_megajy + dust_q_megajy + sync_q_megajy
         flux_u = cmb_u_megajy + dust_u_megajy + sync_u_megajy

;;          flux_i = cmb_t_megajy
;;          flux_q = cmb_q_megajy
;;          flux_u = cmb_u_megajy


         ;; Account for non linearity in measurement
         for ieps = 0, neps-1 do begin 
            m = flux_i + cos2alpha*flux_q + sin2alpha*flux_u
            m_eps = m + epsilon[ieps]*m^2
            
            ;; Convert back to thermodynamic temperature
            convert_megajy_millik, lambda_microns, m_eps, m_mk, /cmb
            m_microK = m_mk*1d3
            
            atd[*,0,ieps] +=  m_microK
            atd[*,1,ieps] +=  m_microK*cos2alpha
            atd[*,2,ieps] +=  m_microK*sin2alpha
            atd[*,3,ieps] += rj2k * (nu[inu]/nu0)^beta_dust           * m_microK
            atd[*,4,ieps] += rj2k * (nu[inu]/nu0)^beta_dust*cos2alpha * m_microK
            atd[*,5,ieps] += rj2k * (nu[inu]/nu0)^beta_dust*sin2alpha * m_microK
            atd[*,6,ieps] += rj2k * (nu[inu]/nu0)^beta_sync           * m_microK
            atd[*,7,ieps] += rj2k * (nu[inu]/nu0)^beta_sync*cos2alpha * m_microK
            atd[*,8,ieps] += rj2k * (nu[inu]/nu0)^beta_sync*sin2alpha * m_microK
         endfor
      endfor
   endfor

   ;; Perform component separation
   for ieps=0, neps-1 do begin 

      s[*,*,ieps] = atam1##atd[*,*,ieps]
      s1 = s[*,*,ieps]
      ;; CMB output is directly is microK_CMB
      out_cmb_t[*,*,ieps]  = reform( s1[*,0], nx, ny)
      out_cmb_q[*,*,ieps]  = reform( s1[*,1], nx, ny)
      out_cmb_u[*,*,ieps]  = reform( s1[*,2], nx, ny)

;;    wind, 1, 1, /free, /large
;;    my_multiplot, 3, 4, pp, pp1, /rev, /full, /dry
;;    dp = {noerase:1, charsize:0.6, charbar:0.6}
;;    imview, nu0_dust_t, position=pp1[0,*], title='Input dust T', dp=dp
;;    imview, nu0_dust_q, position=pp1[1,*], title='Input dust Q', dp=dp
;;    imview, nu0_dust_u, position=pp1[2,*], title='Input dust U', dp=dp
;; 
;;    imview, nu0_sync_t, position=pp1[3,*], title='Input sync T', dp=dp
;;    imview, nu0_sync_q, position=pp1[4,*], title='Input sync Q', dp=dp
;;    imview, nu0_sync_u, position=pp1[5,*], title='Input sync U', dp=dp
;; 
;;    imview, cmb_t, position=pp1[6,*], title='input CMB T', dp=dp
;;    imview, cmb_q, position=pp1[7,*], title='input CMB Q', dp=dp
;;    imview, cmb_u, position=pp1[8,*], title='input CMB U', dp=dp
;; 
;;    imview, out_cmb_t, position=pp1[9,*], title='output T CMB', dp=dp
;;    imview, out_cmb_q, position=pp1[10,*], title='output Q CMB', dp=dp
;;    imview, out_cmb_u, position=pp1[11,*], title='output U CMB', dp=dp
;;

;;    wind, 1, 1, /free, /large
;;    my_multiplot, 3, 3, pp, pp1, /rev
;;    dp = {noerase:1, charsize:0.6, charbar:0.6}
;;    imview, cmb_t, position=pp1[0,*], title='input CMB T', dp=dp
;;    imview, cmb_q, position=pp1[1,*], title='input CMB Q', dp=dp
;;    imview, cmb_u, position=pp1[2,*], title='input CMB U', dp=dp
;; 
;;    imview, out_cmb_t, position=pp1[3,*], title='output T CMB', dp=dp
;;    imview, out_cmb_q, position=pp1[4,*], title='output Q CMB', dp=dp
;;    imview, out_cmb_u, position=pp1[5,*], title='output U CMB', dp=dp
;;
;;    plot, cmb_t, out_cmb_t, psym=1, position=pp1[6,*], /noerase
;;    fit = linfit( cmb_t, out_cmb_t)
;;    oplot, minmax(cmb_t), fit[0] + fit[1]*minmax(cmb_t)
;;    legendastro, ['T', strtrim(fit,2)]
;;
;;    plot, cmb_q, out_cmb_q, psym=1, position=pp1[7,*], /noerase
;;    fit = linfit( cmb_q, out_cmb_q)
;;    oplot, minmax(cmb_q), fit[0] + fit[1]*minmax(cmb_q)
;;    legendastro, ['Q', strtrim(fit,2)]
;;
;;    plot, cmb_u, out_cmb_u, psym=1, position=pp1[8,*], /noerase
;;    fit = linfit( cmb_u, out_cmb_u)
;;    oplot, minmax(cmb_u), fit[0] + fit[1]*minmax(cmb_u)
;;    legendastro, ['U', strtrim(fit,2)]
;;
;;stop

      qu2eb, out_cmb_q[*,*,ieps], out_cmb_u[*,*,ieps], res_arcmin, out_cmb_e, out_cmb_b
      ipoker, out_cmb_t[*,*,ieps], res_arcmin, k, pk_cmb_t, /bypass, /rem, delta_l_over_l=delta_l_over_l
      ipoker, out_cmb_e, res_arcmin, k, pk_cmb_e, /bypass, /rem, delta_l_over_l=delta_l_over_l
      ipoker, out_cmb_b, res_arcmin, k, pk_cmb_b, /bypass, /rem, delta_l_over_l=delta_l_over_l
      ipoker, out_cmb_t[*,*,ieps], res_arcmin, map1=out_cmb_e, k, pk_cmb_te, /bypass, /rem, delta_l_over_l=delta_l_over_l
      
      clt_out_res[ imc,*,ieps] = k*(k+1)/(2*!dpi)*pk_cmb_t
      cle_out_res[ imc,*,ieps] = k*(k+1)/(2*!dpi)*pk_cmb_e
      clb_out_res[ imc,*,ieps] = k*(k+1)/(2*!dpi)*pk_cmb_b
      clte_out_res[imc,*,ieps] = k*(k+1)/(2*!dpi)*pk_cmb_te
   endfor
endfor

clt_out_avg        = dblarr(nk,neps)
cle_out_avg        = dblarr(nk,neps)
clb_out_avg        = dblarr(nk,neps)
clte_out_avg       = dblarr(nk,neps)
sigma_clt_out_avg  = dblarr(nk,neps)
sigma_cle_out_avg  = dblarr(nk,neps)
sigma_clb_out_avg  = dblarr(nk,neps)
sigma_clte_out_avg = dblarr(nk,neps)
frac               = dblarr(neps)

l_ref = 100
k_ref = 100
wk    = (where( abs(k-k_ref) eq min( abs(k-k_ref))))[0]
w     = (where( abs(l-l_ref) eq min( abs(l-l_ref))))[0]

for ieps=0,neps-1 do begin
   mc_reduce, clt_out_res[*,*,ieps], clt_out_avg_res, sigma_clt_out_avg_res
   mc_reduce, cle_out_res[*,*,ieps], cle_out_avg_res, sigma_cle_out_avg_res
   mc_reduce, clb_out_res[*,*,ieps], clb_out_avg_res, sigma_clb_out_avg_res
   mc_reduce, clte_out_res[*,*,ieps], clte_out_avg_res, sigma_clte_out_avg_res

   clt_out_avg[*,ieps]  = clt_out_avg_res
   cle_out_avg[*,ieps]  = cle_out_avg_res
   clb_out_avg[*,ieps]  = clb_out_avg_res
   clte_out_avg[*,ieps] = clte_out_avg_res
   sigma_clt_out_avg[*,ieps] = sigma_clt_out_avg_res
   sigma_cle_out_avg[*,ieps] = sigma_cle_out_avg_res
   sigma_clb_out_avg[*,ieps] = sigma_clb_out_avg_res
   sigma_clte_out_avg[*,ieps] = sigma_clte_out_avg_res

   frac[ieps] = clb_out_avg_res[wk]/clb[w]
endfor



psym = 8
syms = 0.5
;fmt = '(F4.2)'
fmt = '(F5.3)'
png = 0
ps  = 0


eps_min = 5d-3
eps_max = 0.20
neps    = 1000
epsilon1 = dindgen(neps)*(eps_max-eps_min)/(neps-1) + eps_min


a = poly_fit(epsilon, frac, 2)

if ps eq 0 then wind, 1, 1, /free
outplot, file='epsilon_power_law', ps=ps, png=png, thick=thick
plot , epsilon, frac, psym=8, xtitle='!7e!3', ytitle='!12C!dl!n!u!3NL!n / !12C!dl!n!u!3in!n', charsize=1.2
oplot, epsilon1, a[0] + a[1]*epsilon1 + a[2]*epsilon1^2, col=250, psym=4
oplot, [0,10], [1,1], linestyle=2
outplot, /close, /verb

wfit = where( (a[0] + a[1]*epsilon1 + a[2]*epsilon1^2 lt 1.1) and (a[0] + a[1]*epsilon1 + a[2]*epsilon1^2 gt 0.99), nwfit)
print, epsilon1[wfit]

png = 0
ps  = 0 

if ps eq 0 then wind, 1, 1, /free, /large
outplot, file='cmb_power_spectra_nl', ps=ps, png=png, thick=thick
yra = [1d-6, 1d4]
yra = [1d-10,1d4]
xra = [1, max(k)*2]
plot_oo, k, clt_out_avg, xra=xra, /xs, yra=yra, /ys, $
         xtitle='Multipole !12l!3', ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n', $
         psym=psym, syms=syms
oploterror, k, clt_out_avg, sigma_clt_out_avg, psym=psym, syms=syms
oploterror, k, cle_out_avg, sigma_cle_out_avg, psym=psym, syms=syms, col=col_e, errcol=col_e
for ieps=0, neps-1 do begin 
   oploterror, k, clb_out_avg[*,ieps], sigma_clb_out_avg[*,ieps], psym=psym, syms=syms, col=col_b, errcol=col_b
   xyouts, 15, min(clb_out_avg[*,ieps]), '!7e!3 ='+string(epsilon[ieps],form=fmt), charsize=0.8
endfor
oploterror, k, abs(clte_out_avg), sigma_clte_out_avg, psym=psym, syms=syms, col=col_te, errcol=col_te
oplot, l, l*(l+1)/(2*!dpi)*cmb_clt
oplot, l, l*(l+1)/(2*!dpi)*cmb_cle, col=col_e
oplot, l, l*(l+1)/(2*!dpi)*abs(cmb_clte), col=col_te
oplot, l, l*(l+1)/(2*!dpi)*cmb_clb_in, col=col_b
legendastro, ['T', 'E', 'B', 'TE'], col=[0, col_e, col_b, col_te], line=0
xyouts, 2, 5d-5, 'BB (r=0.001)'
outplot, /close, /verb


end
