
;; Simple component separation to assess impact of non linearity

;; Assumed non-linearity parameter
epsilon = 2 ; 0.d0 ; 1d-4 ;0.d0 ;1d-3 ;0.5d-1 ;1d-6 ;0.d0

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
nmc = 3 ; 20 ; 20 ;10 ; 50 ;30

;;============================================
n_nu  = n_elements(nu)
alpha = alpha_det_deg*!dtor
ndet  = n_elements(alpha)

;; Quick check
readcol, '$SK_DIR/Cl/cmb_totCls_r0.001.dat', l, clt, cle, clb, clte, $
         format='D,D,D,D', comment='#', /silent
;; prefactor included in the CAMB's output files
fl = l*(l+1)/(2*!dpi)
cmb_clt  = clt/fl
cmb_cle  = cle/fl
cmb_clb_in = clb/fl
cmb_clte = clte/fl


col_e = 70
col_te = 150
col_b = 250

;; CMB power spectra, from CAMB 
;; ps=0
;; png=0

;; if ps eq 0 then wind, 1, 1, /free, /large
;; outplot, file='cmb_power_spectra', ps=ps, png=png, thick=thick
;; plot, l, l*(l+1)/(2*!dpi)*cmb_clt, /xs, /ys, /xlog, /ylog, yra=[1d-10,1d5], $
;;          xtitle='Multipole !12l!3', ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n'
;; oplot, l, l*(l+1)/(2*!dpi)*cmb_cle, col=col_e
;; oplot, l, l*(l+1)/(2*!dpi)*abs(cmb_clte), col=col_te
;; oplot, l, l*(l+1)/(2*!dpi)*cmb_clb_in, col=col_b
;; legendastro, ['TT', 'EE', 'BB', 'TE'], col=[0, col_e, col_b, col_te], line=0
;; xyouts, 2, 5d-5, 'BB (r=0.001)'
;; outplot, /close, /verb

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
nu0 = 100 ; to be simple in the manuscript
rj2k_nu0 = rj2thermo(nu0)

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
thick = 1 ; 2
col_fit = 250
xra_dust = [1,100]
xra_sync = [1,50]
yra_dust = [1d-4, 1d4]
yra_sync = [1d-6, 1d2]
xtitle='Multipole !12l!3'
ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n'

ps=0
png=0
fmt = '(F4.1)'

if ps eq 0 then wind, 1, 1, /free, /large
outplot, file='dust_sync_power_spectra', ps=ps, png=png, thick=thick
my_multiplot, 4, 2, pp, pp1, /rev, gap_x=0.07

plot_oo, l_out, l_out*(l_out+1)/(2*!dpi)*abs(clt_dust), xra=xra_dust, yra=yra_dust, position=pp1[0,*], psym=psym, syms=syms, $
         xtitle=xtitle, ytitle=ytitle
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_t[0])*l_out^dust_fit_t[1], col=col_fit, thick=thick
legendastro, ['Dust', 'TT']
legendastro, ['fit coeff :',$
              'A='+string(exp(dust_fit_t[0]),form=fmt),$
              '!7a!3='+string(dust_fit_t[1],form=fmt)], /bottom, textcol=250, /left

plot_oo,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(cle_dust), psym=psym, syms=syms, $
        xra=xra_dust, yra=yra_dust, position=pp1[1,*], /noerase, xtitle=xtitle, ytitle=ytitle
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_e[0])*l_out^dust_fit_e[1], col=col_fit, thick=thick
legendastro, ['Dust', 'EE']
legendastro, ['fit coeff :',$
              'A='+string(exp(dust_fit_e[0]),form=fmt),$
              '!7a!3='+string(dust_fit_e[1],form=fmt)], /bottom, textcol=250, /left

plot_oo,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clb_dust), psym=psym, syms=syms, $
        xra=xra_dust, yra=yra_dust, position=pp1[2,*], /noerase, xtitle=xtitle, ytitle=ytitle
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_b[0])*l_out^dust_fit_b[1], col=col_fit, thick=thick
legendastro, ['Dust', 'BB']
legendastro, ['fit coeff :',$
              'A='+string(exp(dust_fit_b[0]),form=fmt),$
              '!7a!3='+string(dust_fit_b[1],form=fmt)], /bottom, textcol=250, /left

plot_oo,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clte_dust), psym=psym, syms=syms,$
        xra=xra_dust, yra=yra_dust, position=pp1[3,*], /noerase, xtitle=xtitle, ytitle=ytitle
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_te[0])*l_out^dust_fit_te[1], col=col_fit, thick=thick
legendastro, ['Dust', 'TE']
legendastro, ['fit coeff :',$
              'A='+string(exp(dust_fit_te[0]),form=fmt),$
              '!7a!3='+string(dust_fit_te[1],form=fmt)], /bottom, textcol=250, /left


plot_oo, l_out, l_out*(l_out+1)/(2*!dpi)*abs(clt_sync), xra=xra_sync, yra=yra_sync, position=pp1[4,*],$
         psym=psym, syms=syms, /noerase, xtitle=xtitle, ytitle=ytitle
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_t[0])*l_out^sync_fit_t[1], col=col_fit, thick=thick
legendastro, ['Synchrotron', 'TT']
legendastro, ['fit coeff :',$
              'A='+string(exp(sync_fit_t[0]),form=fmt),$
              '!7a!3='+string(sync_fit_t[1],form=fmt)], /bottom, textcol=250, /left

plot_oo,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(cle_sync), psym=psym, syms=syms, $
        xra=xra_sync, yra=yra_sync, position=pp1[5,*], /noerase, xtitle=xtitle, ytitle=ytitle
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_e[0])*l_out^sync_fit_e[1], col=col_fit, thick=thick
legendastro, ['Synchrotron', 'EE']
legendastro, ['fit coeff :',$
              'A='+string(exp(sync_fit_e[0]),form=fmt),$
              '!7a!3='+string(sync_fit_e[1],form=fmt)], /bottom, textcol=250, /left


plot_oo,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clb_sync), psym=psym, syms=syms, $
        xra=xra_sync, yra=yra_sync, position=pp1[6,*], /noerase, xtitle=xtitle, ytitle=ytitle
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_b[0])*l_out^sync_fit_b[1], col=col_fit, thick=thick
legendastro, ['Synchrotron', 'BB']
legendastro, ['fit coeff :',$
              'A='+string(exp(sync_fit_b[0]),form=fmt),$
              '!7a!3='+string(sync_fit_b[1],form=fmt)], /bottom, textcol=250, /left


plot_oo,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clte_sync), psym=psym, syms=syms,$
        xra=xra_sync, yra=yra_sync, position=pp1[7,*], /noerase, xtitle=xtitle, ytitle=ytitle
oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_te[0])*l_out^sync_fit_te[1], col=col_fit, thick=thick
legendastro, ['Synchrotron', 'TE']
legendastro, ["A="+string(exp(sync_fit_te[0]),form=fmt), $
              "!7a!3="+string(sync_fit_te[1],form=fmt)], $
             /bottom, textcol=250, /left, col=250

outplot, /close, /verb

;; wind, 1, 1, /free, /large
;; my_multiplot, 2, 1, pp, pp1, /rev
;; plot_oo, l_out, l_out*(l_out+1)/(2*!dpi)*abs(clt_dust), xra=xra, yra=[1d-4, 1d4], position=pp1[0,*], psym=psym, syms=syms
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(cle_dust), col=col_e, psym=psym, syms=syms
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clb_dust), col=col_b, psym=psym, syms=syms
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clte_dust), col=col_te, psym=psym, syms=syms
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_t[0])*l_out^dust_fit_t[1]
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_e[0])*l_out^dust_fit_e[1], col=col_e
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_b[0])*l_out^dust_fit_b[1], col=col_b
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(dust_fit_te[0])*l_out^dust_fit_te[1], col=col_te
;; legendastro, 'Dust'

;; plot_oo, l_out, l_out*(l_out+1)/(2*!dpi)*abs(clt_sync), xra=xra, yra=[1d-10, 1d2], position=pp1[1,*],$
;;          /noerase, psym=psym, syms=syms
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(cle_sync), col=col_e, psym=psym, syms=syms
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clb_sync), col=col_b, psym=psym, syms=syms
;; oplot,   l_out, l_out*(l_out+1)/(2*!dpi)*abs(clte_sync), col=col_te, psym=psym, syms=syms
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_t[0])*l_out^sync_fit_t[1]
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_e[0])*l_out^sync_fit_e[1], col=col_e
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_b[0])*l_out^sync_fit_b[1], col=col_b
;; oplot, l_out, l_out*(l_out+1)/(2*!dpi)*exp(sync_fit_te[0])*l_out^sync_fit_te[1], col=col_te
;; legendastro, 'Sync'

;; dust at nu0
dust_clt  = exp(dust_fit_t[ 0])*l^dust_fit_t[ 1]
dust_cle  = exp(dust_fit_e[ 0])*l^dust_fit_e[ 1]
dust_clb  = exp(dust_fit_b[ 0])*l^dust_fit_b[ 1]
dust_clte = exp(dust_fit_te[0])*l^dust_fit_te[1]

print, "Dust (TT) (uK_RJ): "+strtrim(exp(dust_fit_t[ 0]),2)+" x ell^"+strtrim(dust_fit_t[1],2)
print, "Dust (EE) (uK_RJ): "+strtrim(exp(dust_fit_e[ 0]),2)+" x ell^"+strtrim(dust_fit_e[1],2)
print, "Dust (TE) (uK_RJ): "+strtrim(exp(dust_fit_te[ 0]),2)+" x ell^"+strtrim(dust_fit_te[1],2)
print, "Dust (BB) (uK_RJ): "+strtrim(exp(dust_fit_b[ 0]),2)+" x ell^"+strtrim(dust_fit_b[1],2)


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

fl = l*(l+1)/(2*!dpi)
wind, 1, 1, /free, /large
plot, l, fl*cmb_clt, xra=[1, max(l)], /xlog, /xs, yra=[1d-4, 1d4], /ylog
oplot, l, fl*cmb_cle, col=70
oplot, l, fl*abs(cmb_clte), col=150
oplot, l, fl*dust_cle, col=200
oplot, l, fl*dust_clb, col=200, line=2
oplot, l, fl*sync_cle, col=100
oplot, l, fl*sync_clb, col=100, line=2

;; CMB, dust, synchrotron maps
ps=0
png=1
xymaps, nx, nx, -nx/2*res_arcmin/60., -nx/2*res_arcmin/60., res_arcmin/60., xmap, ymap
if ps eq 0 then wind, 1, 1, /free, xs=1100, ys=800
outplot, file='cmb_dust_sync_maps', ps=ps, png=png, thick=thick
dp = {noerase:1, charsize:0.6, charbar:0.6, xmap:xmap, ymap:ymap}
delvarx, unitsbar, unitsbar1
;unitsbar='[!7l!3K!dcmb!n]'
;unitsbar1 = '[!7l!3K!dRJ!n]'
my_multiplot, 3, 3, pp, pp1, /rev, gap_x=0.05, gap_y=0.05, $
              xmin=0.05, xmax=0.98, ymin=0.02, ymax=0.98, xmargin=0.001, ymargin=0.001
imview, cmb_t, dp=dp, position=pp1[0,*], title='CMB T', unitsbar=unitsbar
imview, cmb_q, dp=dp, position=pp1[1,*], title='CMB Q', unitsbar=unitsbar
imview, cmb_u, dp=dp, position=pp1[2,*], title='CMB U', unitsbar=unitsbar

imview, dust_t * rj2k_nu0, dp=dp, position=pp1[3,*], title='DUST T', unitsbar=unitsbar1
imview, dust_q * rj2k_nu0, dp=dp, position=pp1[4,*], title='DUST Q', unitsbar=unitsbar1
imview, dust_u * rj2k_nu0, dp=dp, position=pp1[5,*], title='DUST U', unitsbar=unitsbar1
                                                                           
imview, sync_t * rj2k_nu0, dp=dp, position=pp1[6,*], title='SYNC T', unitsbar=unitsbar1
imview, sync_q * rj2k_nu0, dp=dp, position=pp1[7,*], title='SYNC Q', unitsbar=unitsbar1
imview, sync_u * rj2k_nu0, dp=dp, position=pp1[8,*], title='SYNC U', unitsbar=unitsbar1
outplot, /close, /verb

;stop



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

;;=================================== Main loop =======================
nk = n_elements(k)
nu = double([70, 100, 150, 217, 353])
nu = double([70, 150, 217])
epsilon_list = [1.d-1, 1.d-2, 1.d-3, 1.d-4, 1d-5]
;epsilon_list = [0.d0]
n_epsilon = n_elements(epsilon_list)
n_nu = n_elements(nu)
clt_out_raw  = dblarr(nmc,n_nu,n_epsilon,nk)
cle_out_raw  = dblarr(nmc,n_nu,n_epsilon,nk)
clb_out_raw  = dblarr(nmc,n_nu,n_epsilon,nk)
clte_out_raw = dblarr(nmc,n_nu,n_epsilon,nk)
cltb_out_raw = dblarr(nmc,n_nu,n_epsilon,nk)
cleb_out_raw = dblarr(nmc,n_nu,n_epsilon,nk)

for imc=0, nmc-1 do begin
   print, strtrim(imc,2)+"/"+strtrim(nmc-1,2)
   ;; cmb spectra are in microK_CMB^2
   cls2mapsiqu, l, cmb_clt, cmb_cle, cmb_clb, cmb_clte, nx, res_arcmin/60., $
                cmb_t, cmb_q, cmb_u
   ;; dust and sync spectra are in microK_RJ^2
   cls2mapsiqu, l, dust_clt, dust_cle, dust_clb, dust_clte, nx, res_arcmin/60., $
                nu0_dust_t, nu0_dust_q, nu0_dust_u
   cls2mapsiqu, l, sync_clt, sync_cle, sync_clb, sync_clte, nx, res_arcmin/60., $
                nu0_sync_t, nu0_sync_q, nu0_sync_u, cu_t, cu_e, cu_b, cu_te

   for inu=0, n_nu-1 do begin
      rj2k = rj2thermo( nu[inu])
      dust_t = rj2k * (nu[inu]/nu0)^beta_dust * nu0_dust_t
      dust_q = rj2k * (nu[inu]/nu0)^beta_dust * nu0_dust_q
      dust_u = rj2k * (nu[inu]/nu0)^beta_dust * nu0_dust_u
      
      sync_t = rj2k * (nu[inu]/nu0)^beta_sync * nu0_sync_t
      sync_q = rj2k * (nu[inu]/nu0)^beta_sync * nu0_sync_q
      sync_u = rj2k * (nu[inu]/nu0)^beta_sync * nu0_sync_u
      
      ;; Effective Stokes parameters
      I = cmb_t + dust_t + sync_t
      Q = cmb_q + dust_q + sync_q
      U = cmb_u + dust_u + sync_u

      for ieps=0, n_epsilon-1 do begin
;;          I_eff = I + epsilon_list[ieps] * I^2
;;          Q_eff = Q + epsilon_list[ieps] * 2.d*I*Q
;;          U_eff = U + epsilon_list[ieps] * 2.d*I*U
         I_eff = I + epsilon_list[ieps] * I^2
         Q_eff = epsilon_list[ieps] * 2.d*I*Q
         U_eff = epsilon_list[ieps] * 2.d*I*U

         qu2eb, Q_eff, U_eff, res_arcmin, E_eff, B_eff
         ipoker, I_eff, res_arcmin, k, pk_cmb_t, /bypass, /rem, delta_l_over_l=delta_l_over_l
         ipoker, E_eff, res_arcmin, k, pk_cmb_e, /bypass, /rem, delta_l_over_l=delta_l_over_l
         ipoker, B_eff, res_arcmin, k, pk_cmb_b, /bypass, /rem, delta_l_over_l=delta_l_over_l
         ipoker, I_eff, res_arcmin, map1=E_eff, k, pk_cmb_te, /bypass, /rem, delta_l_over_l=delta_l_over_l
         ipoker, I_eff, res_arcmin, map1=B_eff, k, pk_cmb_tb, /bypass, /rem, delta_l_over_l=delta_l_over_l
         ipoker, E_eff, res_arcmin, map1=B_eff, k, pk_cmb_eb, /bypass, /rem, delta_l_over_l=delta_l_over_l

         clt_out_raw[imc,inu,ieps,*] = pk_cmb_t
         cle_out_raw[imc,inu,ieps,*] = pk_cmb_e
         clb_out_raw[imc,inu,ieps,*] = pk_cmb_b
         clte_out_raw[imc,inu,ieps,*] = pk_cmb_te
         cltb_out_raw[imc,inu,ieps,*] = pk_cmb_tb
         cleb_out_raw[imc,inu,ieps,*] = pk_cmb_eb
         
      endfor
   endfor
endfor


wind, 1, 1, /free, /large
my_multiplot, 1, 1, ntot=n_nu, pp, pp1, /rev
yra = [1d-7,1d5]
xra = [1, max(k)*2]
thick=1
epsilon_color = [40, 80, 150, 200, 250]
for inu=0, n_nu-1 do begin
   plot_oo, xra, yra, /nodata, xra=xra, /xs, yra=yra, /ys, $
            xtitle='Multipole l', ytitle='l(l+1)C!dl!n/2!7p!3 !7l!3K!u2!n', $
            psym=psym, syms=syms, position=pp1[inu,*], /noerase
   legendastro, strtrim( long(nu[inu]),2)+" GHz", /right
   for ieps=0, n_epsilon-1 do begin
      mc_reduce, reform( clt_out_raw[*,inu,ieps,*], nmc, nk), clt_avg, sigma_clt_avg
      mc_reduce, reform( cle_out_raw[*,inu,ieps,*], nmc, nk), cle_avg, sigma_cle_avg
      mc_reduce, reform( clb_out_raw[*,inu,ieps,*], nmc, nk), clb_avg, sigma_clb_avg
      mc_reduce, reform( clte_out_raw[*,inu,ieps,*], nmc, nk), clte_avg, sigma_clte_avg
      mc_reduce, reform( cltb_out_raw[*,inu,ieps,*], nmc, nk), cltb_avg, sigma_cltb_avg
      mc_reduce, reform( cleb_out_raw[*,inu,ieps,*], nmc, nk), cleb_avg, sigma_cleb_avg

;;      oploterror, k, k*(k+1)/(2*!dpi)*clt_avg,  k*(k+1)/(2*!dpi)*sigma_clt_avg, psym=8, syms=0.5
;      oploterror, k, k*(k+1)/(2*!dpi)*cle_avg,  k*(k+1)/(2*!dpi)*sigma_cle_avg, psym=8, syms=0.5
      oplot, k, k*(k+1)/(2*!dpi)*clb_avg, psym=-8, syms=0.5, col=epsilon_color[ieps]
;      oploterror, k, k*(k+1)/(2*!dpi)*clte_avg, k*(k+1)/(2*!dpi)*sigma_clte_avg, psym=8, syms=0.5
;      oploterror, k, k*(k+1)/(2*!dpi)*cltb_avg, k*(k+1)/(2*!dpi)*sigma_cltb_avg, psym=8, syms=0.5
;      oploterror, k, k*(k+1)/(2*!dpi)*cleb_avg, k*(k+1)/(2*!dpi)*sigma_cleb_avg, psym=8, syms=0.5
      legendastro, '!7e!3 = '+strtrim( string(epsilon_list,form='(F6.4)'),2), /bottom, $
                   col=epsilon_color
   endfor
   oplot, l, l*(l+1)/(2*!dpi)*cmb_clt
   oplot, l, l*(l+1)/(2*!dpi)*cmb_cle
   oplot, l, l*(l+1)/(2*!dpi)*abs(cmb_clte)
   oplot, l, l*(l+1)/(2*!dpi)*cmb_clb_in, line=2
   legendastro, ['T', 'E', 'B', 'TE'], line=[0,0,2,0]
endfor

;; Plots
;; Pour le manuscrit
png = 1
epsilon_limit = dblarr(n_nu)
wind, 1, 1, /free, /large
for inu=0, n_nu-1 do begin
   outplot, file='cls_raw_epsilon_'+strtrim(long(nu[inu]),2)+"GHz", png=png
   plot_oo, xra, yra, /nodata, xra=xra, /xs, yra=yra, /ys, $
            xtitle='Multipole l', ytitle='l(l+1)C!dl!n/2!7p!3 !7l!3K!u2!n', $
            psym=psym, syms=syms
   legendastro, strtrim( long(nu[inu]),2)+" GHz", /right
   for ieps=0, n_epsilon-1 do begin
      mc_reduce, reform( clb_out_raw[*,inu,ieps,*], nmc, nk), clb_avg, sigma_clb_avg
      oplot, k, k*(k+1)/(2*!dpi)*clb_avg, psym=-8, syms=0.5, col=epsilon_color[ieps]
      if ieps eq 0 then clb_avg_100 = interpol( k*(k+1)/(2*!dpi)*clb_avg, k, 100.)
   endfor
   oplot, l, l*(l+1)/(2*!dpi)*cmb_clt
   oplot, l, l*(l+1)/(2*!dpi)*cmb_cle
   oplot, l, l*(l+1)/(2*!dpi)*abs(cmb_clte)
   oplot, l, l*(l+1)/(2*!dpi)*cmb_clb_in, line=2
   legendastro, ['T', 'E', 'B', 'TE'], line=[0,0,2,0]
   legendastro, '!7e!3 = 10!u'+strtrim( long(alog10(epsilon_list)),2)+"!n", $
                /bottom, col=epsilon_color
   outplot, /close, /verb

   ;; request a 1/10th of the target B mode at l=100 (r=0.001)
   ieps = 0
   w100 = (where( l eq 100))[0]
   epsilon_limit[inu] = sqrt( (l*(l+1)/(2*!dpi)*cmb_clb_in)[w100]/10./clb_avg_100 * epsilon_list[ieps]^2)
endfor

print, ""
for inu=0, n_nu-1 do print, ";; Espilon limit ("+strtrim( long(nu[inu]),2)+" GHz): "+strtrim(epsilon_limit[inu],2)
;; Espilon limit (70 GHz): 2.2458870e-05
;; Espilon limit (150 GHz): 1.9674821e-05
;; Espilon limit (217 GHz): 9.7967775e-06


end
