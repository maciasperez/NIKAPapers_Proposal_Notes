
;; Simple component separation to assess impact of non linearity

;; Assumed non-linearity parameter
epsilon = 0.5d-1

;pro a, epsilon

;; Outplot format
ps  = 0
png = 1

;; Final map resolution
nside = 256 ; 1024 ; 512
noplot = 0

;; Foreground spectral parameters
beta_dust = 1.5
beta_sync = -3

;; Experiment bands
nu = [70.d0, 100.d0, 143.d0, 217.d0, 353.d0];, 545.d0]

;; Take a n-detector system with perfect angular coverage per band
alpha_det_deg = [0.d0, 60.d0, 120.d0]

;; MC realizations
nmc = 30

;; Latitude mask
latitude_cut = 20 ; 30 ; deg

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

   ata = dblarr(9,9)
   atd = dblarr(npix,9)
   
   for inu=0, n_nu-1 do begin
      rj2k = rj2thermo(nu[inu])
      
      for ialpha=0, n_alpha-1 do begin
         cos2alpha = cos(2.*alpha[ialpha])
         sin2alpha = sin(2.*alpha[ialpha])
         
         ;; Fill ata
         ata[0,0] += 1.d0
         ata[1,0] += cos2alpha
         ata[2,0] += sin2alpha
         ata[3,0] += rj2k * (nu[inu]/dust_i_nu0)^beta_dust
         ata[4,0] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha
         ata[5,0] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha
         ata[6,0] += rj2k * (nu[inu]/sync_i_nu0)^beta_sync
         ata[7,0] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha
         ata[8,0] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha

         ata[1,1] += cos2alpha^2
         ata[2,1] += cos2alpha*sin2alpha
         ata[3,1] += rj2k * (nu[inu]/dust_i_nu0)^beta_dust*cos2alpha
         ata[4,1] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha^2
         ata[5,1] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha*cos2alpha
         ata[6,1] += rj2k * (nu[inu]/sync_i_nu0)^beta_sync*cos2alpha
         ata[7,1] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha^2
         ata[8,1] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha*cos2alpha

         ata[2,2] += sin2alpha^2
         ata[3,2] += rj2k * (nu[inu]/dust_i_nu0)^beta_dust*sin2alpha
         ata[4,2] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha*sin2alpha
         ata[5,2] += rj2k * (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha^2
         ata[6,2] += rj2k * (nu[inu]/sync_i_nu0)^beta_sync*sin2alpha
         ata[7,2] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha*sin2alpha
         ata[8,2] += rj2k * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha^2

         ata[3,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^(2.*beta_dust)
         ata[4,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha
         ata[5,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha
         ata[6,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/sync_i_nu0)^beta_sync
         ata[7,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha
         ata[8,3] += rj2k^2 * (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha

         ata[4,4] += rj2k^2 * (nu[inu]/dust_p_nu0)^(2.*beta_dust)*cos2alpha^2
         ata[5,4] += rj2k^2 * (nu[inu]/dust_p_nu0)^(2.*beta_dust)*cos2alpha*sin2alpha
         ata[6,4] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_i_nu0)^beta_sync*cos2alpha
         ata[7,4] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha^2
         ata[8,4] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha*sin2alpha

         ata[5,5] += rj2k^2 * (nu[inu]/dust_p_nu0)^(2.*beta_dust)*sin2alpha^2
         ata[6,5] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_i_nu0)^beta_sync*sin2alpha
         ata[7,5] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha*sin2alpha
         ata[8,5] += rj2k^2 * (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha^2

         ata[6,6] += rj2k^2 * (nu[inu]/sync_i_nu0)^(2.*beta_sync)
         ata[7,6] += rj2k^2 * (nu[inu]/sync_i_nu0)^beta_sync * (nu[inu]/sync_p_nu0)^beta_sync* cos2alpha
         ata[8,6] += rj2k^2 * (nu[inu]/sync_i_nu0)^beta_sync * (nu[inu]/sync_p_nu0)^beta_sync* sin2alpha

         ata[7,7] += rj2k^2 * (nu[inu]/sync_p_nu0)^(2.*beta_sync)*cos2alpha^2
         ata[8,7] += rj2k^2 * (nu[inu]/sync_p_nu0)^(2.*beta_sync)*cos2alpha*sin2alpha

         ata[8,8] += rj2k^2 * (nu[inu]/sync_p_nu0)^(2.*beta_sync)*sin2alpha^2

         ;;------------------------------
         ;; Define A^tD
         
;;         ;; Scale components to the current frequency in microK_RJ then to microK_CMB
;;         i_dust = rj2k*(nu[inu]/dust_i_nu0)^beta_dust*i_dust_ref
;;         q_dust = rj2k*(nu[inu]/dust_p_nu0)^beta_dust*q_dust_ref
;;         u_dust = rj2k*(nu[inu]/dust_p_nu0)^beta_dust*u_dust_ref
;;         i_sync = rj2k*(nu[inu]/sync_i_nu0)^beta_sync*i_sync_ref
;;         q_sync = rj2k*(nu[inu]/sync_p_nu0)^beta_sync*q_sync_ref
;;         u_sync = rj2k*(nu[inu]/sync_p_nu0)^beta_sync*u_sync_ref

         ;; Scale components to the current frequency in microK_RJ
         i_dust = (nu[inu]/dust_i_nu0)^beta_dust*i_dust_ref
         q_dust = (nu[inu]/dust_p_nu0)^beta_dust*q_dust_ref
         u_dust = (nu[inu]/dust_p_nu0)^beta_dust*u_dust_ref
         i_sync = (nu[inu]/sync_i_nu0)^beta_sync*i_sync_ref
         q_sync = (nu[inu]/sync_p_nu0)^beta_sync*q_sync_ref
         u_sync = (nu[inu]/sync_p_nu0)^beta_sync*u_sync_ref

         ;; Convert to Mjy/sr
         lambda_microns = !const.c/(nu[inu]*1d9)*1d6
         convert_millik_megajy, lambda_microns, i_dust/1000., i_dust_megajy, /rj
         convert_millik_megajy, lambda_microns, q_dust/1000., q_dust_megajy, /rj
         convert_millik_megajy, lambda_microns, u_dust/1000., u_dust_megajy, /rj
         convert_millik_megajy, lambda_microns, i_sync/1000., i_sync_megajy, /rj
         convert_millik_megajy, lambda_microns, q_sync/1000., q_sync_megajy, /rj
         convert_millik_megajy, lambda_microns, u_sync/1000., u_sync_megajy, /rj

         convert_millik_megajy, lambda_microns, i_cmb/1000., i_cmb_megajy, /cmb
         convert_millik_megajy, lambda_microns, q_cmb/1000., q_cmb_megajy, /cmb
         convert_millik_megajy, lambda_microns, u_cmb/1000., u_cmb_megajy, /cmb

         ;; Derive measured flux in Mjy
         fi = (i_cmb_megajy + i_dust_megajy + i_sync_megajy)
         fq = (q_cmb_megajy + q_dust_megajy + q_sync_megajy)
         fu = (u_cmb_megajy + u_dust_megajy + u_sync_megajy)

         ;; Account for non linearity
         m = fi + cos2alpha*fq + sin2alpha*fu +$
             epsilon * (fi + cos2alpha*fq + sin2alpha*fu)^2

         ;; Convert back to thermodynamic temperature
         convert_megajy_millik, lambda_microns, m, m_mk, /cmb
         m_microK = m_mk*1d3
         
         atd[*,0] += m_microK
         atd[*,1] += m_microK*cos2alpha
         atd[*,2] += m_microK*sin2alpha
         atd[*,3] += rj2k * m_microK*(nu[inu]/dust_i_nu0)^beta_dust
         atd[*,4] += rj2k * m_microK*(nu[inu]/dust_p_nu0)^beta_dust*cos2alpha
         atd[*,5] += rj2k * m_microK*(nu[inu]/dust_p_nu0)^beta_dust*sin2alpha
         atd[*,6] += rj2k * m_microK*(nu[inu]/sync_i_nu0)^beta_sync
         atd[*,7] += rj2k * m_microK*(nu[inu]/sync_p_nu0)^beta_sync*cos2alpha
         atd[*,8] += rj2k * m_microK*(nu[inu]/sync_p_nu0)^beta_sync*sin2alpha
         
      endfor
   endfor
;; complete ata
   for i=1, 8 do begin
      for j=0, i-1 do begin
         ata[j,i] = ata[i,j]
      endfor
   endfor

   s = invert(ata)##atd
   out_cmb_i  = reform( s[*,0])
   out_cmb_q  = reform( s[*,1])
   out_cmb_u  = reform( s[*,2])
;;    out_dust_i = reform( s[*,3])
;;    out_dust_q = reform( s[*,4])
;;    out_dust_u = reform( s[*,5])
;;    out_sync_i = reform( s[*,6])
;;    out_sync_q = reform( s[*,7])
;;    out_sync_u = reform( s[*,8])
;; 
;;   mollview, cover*(out_cmb_i-i_cmb), title='CMB I out-in, epsilon='+strtrim(epsilon,2)
;; mollview, cover*(out_cmb_q-q_cmb), title='CMB Q out-in, epsilon='+strtrim(epsilon,2)
;;   stop
;; ;; mollview, out_cmb_q-q_cmb, title='CMB Q out-in'
;; ;; mollview, out_dust_i-i_dust_ref, title='Dust (ref nu) out-in'
;; ;; mollview, out_sync_q-q_sync_ref, title='Sync (ref nu) out-in (RJ)'
;; ;; print, "HERE"
;; ;; stop

   ;; Power spectra
   xmap = dblarr(npix,3)
   xmap[*,0] = out_cmb_i
   xmap[*,1] = out_cmb_q
   xmap[*,2] = out_cmb_u
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

;;    ;; One realization
;;    fl = l_out*(l_out+1)/(2.*!dpi)
;;    xra = [1,3*nside]
;;    yra = [1d-8,1d4]
;;    col_ee = 70
;;    col_bb = 250
;;    col_te = 150
;;    wind, 1, 1, /free, /large
;;    plot, xra, yra, xtitle='!12l!3', ytitle='!12l(l+1)C!dl!n/2!7p!3', $
;;          /nodata, /xlog, /ylog, /xs, /ys
;;    legendastro, 'Epsilon '+strtrim(epsilon,2), /bottom
;;    oplot, l, clt
;;    oplot, l, cle
;;    oplot, l, clb
;;    oplot, l, abs(clte)
;;    oplot, l_out, fl*clt_out
;;    oplot, l_out, fl*cle_out, col=col_ee
;;    oplot, l_out, fl*clb_out, col=col_bb
;;    oplot, l_out, fl*abs(clte_out), col=col_te
;;    legendastro, ['T', 'E', 'B', 'TE'], col=[!p.color, col_ee, col_bb, col_te], line=0
;;    stop
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
