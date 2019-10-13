
;; Simple component separation to assess impact of non linearity

;; Outplot format
ps  = 0
png = 0

;; Final map resolution
nside = 256 ;256 ; 1024 ; 512
noplot = 0

;; Foreground spectral parameters
beta_dust = 1.5
beta_sync = -3

;; Experiment bands
nu = [70.d0, 100.d0, 143.d0, 217.d0, 353.d0];, 545.d0]

;; Take a n-detector system with perfect angular coverage per band
alpha_det_deg = [0.d0, 60.d0, 120.d0]

;; MC realizations
nmc = 500 ;30

;; Latitude mask
latitude_cut = 20 ;30 ; deg

;;============================================
n_nu = n_elements(nu)
alpha = alpha_det_deg*!dtor
n_alpha = n_elements(alpha)

;; Restore Planck maps in microK_RJ (except the cmb maps)
noplot = 1
if defined(i_cmb) eq 0 then begin
   get_planck_maps, nside, $
                    i_cmb, q_cmb, u_cmb, $
                    i_dust_ref, q_dust_ref, u_dust_ref, $
                    i_sync_ref, q_sync_ref, u_sync_ref, noplot=noplot
endif


;; no mask (OK)
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

;;-----------------------------------------------
;; Replace COMMANDER CMB Maps by pure simulations
readcol, '$SK_DIR/Cl/cmb_totCls_r0.001.dat', l, clt, cle, clb, clte, $
         format='D,D,D,D', comment='#'
;; cancel BB
;; Need to add monopole and dipole for synfast
cl_in = dblarr(max(l)+1,4)
;; prefactor included in the CAMB's output files
fl = l*(l+1)/(2*!dpi)
cl_in[2:*,0] = clt/fl
cl_in[2:*,1] = cle/fl
cl_in[2:*,2] = 0.d0 ; clb/fl
cl_in[2:*,3] = clte/fl
input_cl_file = 'mycl.fits'
cl2fits, cl_in, input_cl_file
npix = nside2npix( nside)
cl_in = mrdfits( input_cl_file, 1, cl_header)
ncl_in = n_elements(cl_in)

dust_p_nu0 = 353.d0
dust_i_nu0 = 545.d0
sync_p_nu0 = 30.d0
sync_i_nu0 = 0.408d0

out_cmb_i  = dblarr(npix)
out_cmb_q  = dblarr(npix)
out_cmb_u  = dblarr(npix)
out_dust_i = dblarr(npix)
out_dust_q = dblarr(npix)
out_dust_u = dblarr(npix)
out_sync_i = dblarr(npix)
out_sync_q = dblarr(npix)
out_sync_u = dblarr(npix)

for imc=0, nmc-1 do begin
   isynfast, input_cl_file, cmb_maps, nlmax=3*nside-1, nside=nside
   i_cmb = reform( cmb_maps[*,0])
   q_cmb = reform( cmb_maps[*,1])
   u_cmb = reform( cmb_maps[*,2])

   ;; Power spectra
   xmap = dblarr(npix,3)
   xmap[*,0] = i_cmb
   xmap[*,1] = q_cmb
   xmap[*,2] = u_cmb
   ispice, xmap, cover, l_out, clt_out, cle_out, clb_out, clte_out
   l_out    = l_out[   2:*]
   clt_out  = clt_out[ 2:*]
   cle_out  = cle_out[ 2:*]
   clb_out  = clb_out[ 2:*]
   clte_out = clte_out[2:*]
   
;; Correct for Healpix's pixel window function
   wl = mrdfits( !healpix.path.data+"/pixel_window_n"+string(nside,form='(I4.4)')+".fits",1)
   clt_out  = clt_out[ 2:*]/wl.temperature^2
   cle_out  = cle_out[ 2:*]/wl.polarization^2
   clb_out  = clb_out[ 2:*]/wl.polarization^2
   clte_out = clte_out[2:*]/wl.temperature/wl.polarization

   if imc eq 0 then begin
      clt_out_res  = dblarr(nmc,n_elements(clt_out))
      cle_out_res  = clt_out_res
      clb_out_res  = clt_out_res
      clte_out_res = clt_out_res
   endif
   clt_out_res[imc,*]  = clt_out
   cle_out_res[imc,*]  = cle_out
   clb_out_res[imc,*]  = clb_out
   clte_out_res[imc,*] = clte_out
endfor

;; Plot average results
mc_reduce, clt_out_res, clt_out_avg, sigma_clt_out_avg
mc_reduce, cle_out_res, cle_out_avg, sigma_cle_out_avg
mc_reduce, clb_out_res, clb_out_avg, sigma_clb_out_avg
mc_reduce, clte_out_res, clte_out_avg, sigma_clte_out_avg

;; Derive constraint on epsilon
l_ref = 100
w_out = (where( abs(l_out-l_ref) eq min( abs(l_out-l_ref))))[0]
w     = (where( abs(l-l_ref) eq min( abs(l-l_ref))))[0]
print, "epsilon, clb_out_avg(l="+strtrim(l_ref,2)+")/cl_cmb = ", $
       clb_out_avg[w]/(clb[w]*2*!dpi/(l[w]*(l[w]+1)))

;; Final plot
fl = l_out*(l_out+1)/(2.*!dpi)
xra = [1,3*nside]
yra = [1d-8,1d4]
col_ee = 70
col_bb = 250
col_te = 150
wind, 1, 1, /free, /large
outplot, file='non_linear_fg_residuals_cl_epsilon_'+strtrim(epsilon,2), png=png, ps=ps
plot, xra, yra, xtitle='!12l!3', ytitle='!12l(l+1)C!dl!n/2!7p!3', $
      /nodata, /xlog, /ylog, /xs, /ys
legendastro, 'Epsilon '+strtrim(epsilon,2), /bottom
oplot, l, clt
oplot, l, cle
oplot, l, clb
oplot, l, abs(clte)
oplot, l_out, fl*clt_out_avg
oplot, l_out, fl*cle_out_avg, col=col_ee
oplot, l_out, fl*clb_out_avg, col=col_bb
oplot, l_out, fl*abs(clte_out_avg), col=col_te
legendastro, ['T', 'E', 'B', 'TE'], col=[!p.color, col_ee, col_bb, col_te], line=0
legendastro, 'Nmc = '+strtrim(nmc,2), /right
outplot, /close, /verb


end
