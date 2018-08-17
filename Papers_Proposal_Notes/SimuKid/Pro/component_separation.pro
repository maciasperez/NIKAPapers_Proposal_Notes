
;; Simple component separation to assess impact of non linearity

;; Outplot format
ps  = 0
png = 0

;; Final map resolution
nside = 32 ; 256 ; 1024 ; 512
noplot = 0

;; Foreground spectral parameters
beta_dust = 1.5
beta_sync = -3

;; Experiment bands
nu = [70.d0, 100.d0, 143.d0, 217.d0, 353.d0, 545.d0]

;; Take a n-detector system with perfect angular coverage per band
alpha_det_deg = [0.d0, 60.d0, 120.d0]

;; Assumed non-linearity parameter
epsilon = 0.d0


;;============================================
n_nu = n_elements(nu)
alpha = alpha_det_deg*!dtor
n_alpha = n_elements(alpha)

;; Restore Planck maps and converts them all to microK_CMB
get_planck_maps, nside, $
                 i_cmb, q_cmb, u_cmb, $
                 i_dust_ref, q_dust_ref, u_dust_ref, $
                 i_sync_ref, q_sync_ref, u_sync_ref, noplot=noplot
dust_p_nu0 = 353.d0
dust_i_nu0 = 545.d0
sync_p_nu0 = 30.d0
sync_i_nu0 = 0.408d0

npix = n_elements(i_cmb)
out_cmb_i  = dblarr(npix)
out_cmb_q  = dblarr(npix)
out_cmb_u  = dblarr(npix)
out_dust_i = dblarr(npix)
out_dust_q = dblarr(npix)
out_dust_u = dblarr(npix)
out_sync_i = dblarr(npix)
out_sync_q = dblarr(npix)
out_sync_u = dblarr(npix)

;; Loop over pixels, assuming white noise and "perfect" measures per pix
npix = n_elements(i_cmb)
for ipix=0L, npix-1 do begin

   ;; A^tA matrix and A^td
   ata = dblarr(9,9)
   atd = dblarr(9)
   
   for inu=0, n_nu-1 do begin
      for ialpha=0, n_alpha-1 do begin
         fi = i_cmb[ipix] + $
              (nu[inu]/dust_i_nu0)^beta_dust*i_dust_ref[ipix] + $
              (nu[inu]/sync_i_nu0)^beta_sync*i_sync_ref[ipix]
         fq = q_cmb[ipix] + $
              (nu[inu]/dust_p_nu0)^beta_dust*q_dust_ref[ipix] + $
              (nu[inu]/sync_p_nu0)^beta_sync*q_sync_ref[ipix]
         fu = u_cmb[ipix] + $
              (nu[inu]/dust_p_nu0)^beta_dust*u_dust_ref[ipix] + $
              (nu[inu]/sync_p_nu0)^beta_sync*u_sync_ref[ipix]
              

         cos2alpha = cos(2.*alpha[ialpha])
         sin2alpha = sin(2.*alpha[ialpha])
         
         m = fi + cos2alpha*fq + sin2alpha*fu +$
             epsilon * (fi + cos2alpha*fq + sin2alpha*fu)^2
         
         ata[0,0] += 1.d0
         ata[1,0] += cos2alpha
         ata[2,0] += sin2alpha
         ata[3,0] += (nu[inu]/dust_i_nu0)^beta_dust
         ata[4,0] += (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha
         ata[5,0] += (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha
         ata[6,0] += (nu[inu]/sync_i_nu0)^beta_sync
         ata[7,0] += (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha
         ata[8,0] += (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha

         ata[1,1] += cos2alpha^2
         ata[2,1] += cos2alpha*sin2alpha
         ata[3,1] += (nu[inu]/dust_i_nu0)^beta_dust*cos2alpha
         ata[4,1] += (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha^2
         ata[5,1] += (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha*cos2alpha
         ata[6,1] += (nu[inu]/sync_i_nu0)^beta_sync*cos2alpha
         ata[7,1] += (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha^2
         ata[8,1] += (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha*cos2alpha

         ata[2,2] += sin2alpha^2
         ata[3,2] += (nu[inu]/dust_i_nu0)^beta_dust*sin2alpha
         ata[4,2] += (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha*sin2alpha
         ata[5,2] += (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha^2
         ata[6,2] += (nu[inu]/sync_i_nu0)^beta_sync*sin2alpha
         ata[7,2] += (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha*sin2alpha
         ata[8,2] += (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha^2

         ata[3,3] += (nu[inu]/dust_i_nu0)^(2.*beta_dust)
         ata[4,3] += (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/dust_p_nu0)^beta_dust*cos2alpha
         ata[5,3] += (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/dust_p_nu0)^beta_dust*sin2alpha
         ata[6,3] += (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/sync_i_nu0)^beta_sync
         ata[7,3] += (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha
         ata[8,3] += (nu[inu]/dust_i_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha

         ata[4,4] += (nu[inu]/dust_p_nu0)^(2.*beta_dust)*cos2alpha^2
         ata[5,4] += (nu[inu]/dust_p_nu0)^(2.*beta_dust)*cos2alpha*sin2alpha
         ata[6,4] += (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_i_nu0)^beta_sync*cos2alpha
         ata[7,4] += (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha^2
         ata[8,4] += (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha*sin2alpha

         ata[5,5] += (nu[inu]/dust_p_nu0)^(2.*beta_dust)*sin2alpha^2
         ata[6,5] += (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_i_nu0)^beta_sync*sin2alpha
         ata[7,5] += (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*cos2alpha*sin2alpha
         ata[8,5] += (nu[inu]/dust_p_nu0)^beta_dust * (nu[inu]/sync_p_nu0)^beta_sync*sin2alpha^2

         ata[6,6] += (nu[inu]/sync_i_nu0)^(2.*beta_sync)
         ata[7,6] += (nu[inu]/sync_i_nu0)^beta_sync * (nu[inu]/sync_p_nu0)^beta_sync* cos2alpha
         ata[8,6] += (nu[inu]/sync_i_nu0)^beta_sync * (nu[inu]/sync_p_nu0)^beta_sync* sin2alpha

         ata[7,7] += (nu[inu]/sync_p_nu0)^(2.*beta_sync)*cos2alpha^2
         ata[8,7] += (nu[inu]/sync_p_nu0)^(2.*beta_sync)*cos2alpha*sin2alpha

         ata[8,8] += (nu[inu]/sync_p_nu0)^(2.*beta_sync)*sin2alpha^2


         atd[0] += m
         atd[1] += m*cos2alpha
         atd[2] += m*sin2alpha
         atd[3] += m*(nu[inu]/dust_i_nu0)^beta_dust
         atd[4] += m*(nu[inu]/dust_p_nu0)^beta_dust*cos2alpha
         atd[5] += m*(nu[inu]/dust_p_nu0)^beta_dust*sin2alpha
         atd[6] += m*(nu[inu]/sync_i_nu0)^beta_sync
         atd[7] += m*(nu[inu]/sync_p_nu0)^beta_sync*cos2alpha
         atd[8] += m*(nu[inu]/sync_p_nu0)^beta_sync*sin2alpha

      endfor
   endfor
   
   ;; complete ata
   for i=1, 8 do begin
      for j=0, i-1 do begin
         ata[j,i] = ata[i,j]
;         print, "ata["+strtrim(j,2)+","+strtrim(i,2)+"] = ata["+strtrim(i,2)+","+strtrim(j,2)+"]
      endfor
   endfor
;   print, ata

   s = invert(ata)##atd
   out_cmb_i[ipix]  = s[0]
   out_cmb_q[ipix]  = s[1]
   out_cmb_u[ipix]  = s[2]
   out_dust_i[ipix] = s[3]
   out_dust_q[ipix] = s[4]
   out_dust_u[ipix] = s[5]
   out_sync_i[ipix] = s[6]
   out_sync_q[ipix] = s[7]
   out_sync_u[ipix] = s[8]
   
endfor

mollview, out_cmb_i-i_cmb, title='CMB I out-in'
mollview, out_cmb_q-q_cmb, title='CMB Q out-in'
mollview, out_dust_i-i_dust_ref, title='Dust (ref nu) out-in'
mollview, out_sync_q-q_sync_ref, title='Sync (ref nu) out-in'
print, "HERE"
stop

;; 
;; ;; Cross component maps !
;; ;; amplitude of the foreground residuals
;; xmap = dblarr(npix,3)
;; 
;; ;; Mask out the galactic plane
;; latitude_cut = 30 ; deg
;; cover        = dblarr(npix)
;; ipring       = lindgen(npix)
;; pix2ang_ring, nside, ipring, theta, phi
;; latitude = 90.-theta*!radeg
;; w = where( abs(latitude) ge latitude_cut, nw)
;; cover[w] = 1.d0
;; fsky = float(nw)/npix
;; 
;; ;; Cosmic Variance on BB
;; cosmic_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.curl
;; 
;; ;; Foreground residuals
;; wd, /a
;; ;; Caracterize residual foregrounds in estimated CMB maps by eta
;; eta_pow   = [-2, -3, -4, -5, -6]
;; eta_list   = 10.d0^eta_pow
;; nu_list = [70.d0, 100.d0, 143.d0, 217.d0, 353.d0]
;; 
;; ;eta_list = eta_list[0:1]
;; ;nu_list = [100., 200.]
;; n_eta = n_elements(eta_list)
;; n_nu  = n_elements(nu_list)
;; 
;; epsilon_list = dblarr(n_eta,n_nu)
;; 
;; ;; ianafast, xmap, cl_x, /double, nlmax=3*nside-1, /polar, /ring
;; ;; lx = dindgen(n_elements(cl_x[*,0]))
;; ;; lx = lx[2:*]
;; ;; clxt  = reform( cl_x[2:*,0])
;; ;; clxe  = reform( cl_x[2:*,1])
;; ;; clxb  = reform( cl_x[2:*,2])
;; ;; clxte = reform( cl_x[2:*,3])
;; ;; clxtb = reform( cl_x[2:*,4])
;; ;; clxeb = reform( cl_x[2:*,5])
;; 
;; make_ct, n_eta, mycol
;; fg_res_file = 'fg_res_nl.dat'
;; openw, lu, fg_res_file, /get_lun
;; printf, lu, "#nu (GHz), eta, epsilonT, epsilonE, epsilonB, epsilonTE"
;; for inu=0, n_nu-1 do begin
;; ;   outplot, file='fg_res_spectra_nu'+strtrim( long(nu_list[inu]),2), png=png, ps=ps
;;    plot_init = 0
;; 
;;    ;; Scale foregrounds to the reference frequencies
;;    i_dust = i_dust_ref * (nu_list[inu]/545.d0)^beta_dust
;;    q_dust = q_dust_ref * (nu_list[inu]/353.d0)^beta_dust
;;    u_dust = u_dust_ref * (nu_list[inu]/353.d0)^beta_dust
;; 
;;    i_sync = i_sync_ref * (nu_list[inu]/0.408)^beta_sync
;;    q_sync = q_sync_ref * (nu_list[inu]/30.d0)^beta_sync
;;    u_sync = u_sync_ref * (nu_list[inu]/30.d0)^beta_sync
;; 
;;    ;; Loop on residuals amplitude
;;    for i_eta=0, n_eta-1 do begin
;;       outplot, file='fg_res_spectra_nu'+strtrim( long(nu_list[inu]),2)+$
;;                'eta10'+strtrim(eta_pow[i_eta],2), png=png, ps=ps
;;       eta = eta_list[i_eta]
;;       
;;       xmap[*,0] = cmb_maps[*,0]^2 + $
;;                   2.*cmb_maps[*,0]*eta*(i_dust + i_sync) + $
;;                   eta^2*(i_dust + i_sync)^2
;;       xmap[*,1] = 2.d0*(cmb_maps[*,0]*cmb_maps[*,1] + $
;;                         eta*(q_dust+q_sync)*cmb_maps[*,0] + $
;;                         eta*(i_dust+i_sync)*cmb_maps[*,1])
;;       xmap[*,2] = 2.d0*(cmb_maps[*,0]*cmb_maps[*,2] + $
;;                         eta*(u_dust+u_sync)*cmb_maps[*,0] + $
;;                         eta*(i_dust+i_sync)*cmb_maps[*,2])
;; 
;;       ispice, xmap, cover, lx, clxt, clxe, clxb, clxte
;;       lx    = lx[2:*]
;;       clxt  = clxt[2:*]
;;       clxe  = clxe[2:*]
;;       clxb  = clxb[2:*]
;;       clxte = clxte[2:*]
;;       
;; ;; Pixel window functions
;;       wlpix = mrdfits( !nika.soft_dir+"/idl-libs/Healpix/data/pixel_window_n"+$
;;                        string( nside, format='(I4.4)')+".fits", 1)
;; ;; remove monopole and dipole
;;       wlpix = wlpix[2:ncl_x-1]
;;       
;; ;; Derive spec
;;       w100 = where( l_in eq 100)
;;       wx100 = where( lx eq 100)
;; ;      tol = 0.1
;;       
;;       if plot_init eq 0 then begin
;; ;; All spectra on the same plot
;;          yra = [1d-7, 1d15]
;;          xra = [2, 3000]
;;          thick=2
;;          if ps eq 0 then wind, 1, 1, /free, /large
;;          plot, lx, fl_x*clxt, /xs, /ys, xra=xra, yra=yra, $
;;                xtitle='!12l', $
;;                ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n', $
;;                /ylog, /xlog, /nodata
;;          xyouts, 2.2, (fl_in*cl_in.temperature)[2]*1.3, 'TT'
;;          xyouts, 2.2, (fl_in*cl_in.gradient)[2]*1.3, 'EE'
;;          xyouts, 2.2, (fl_in*cl_in.curl)[2]*1.3, 'BB (r=0.001)'
;;          xyouts, 2.2, (fl_in*abs(cl_in.g_t))[2]*1.3, 'TE'
;; ;         oplot, lx, fl_x*clxt, col=xcol[0]
;;          oplot, l_in, fl_in*cl_in.temperature, col=0
;;          oplot, l_in, fl_in*cl_in.gradient, col=cmb_col
;;          oplot, l_in, fl_in*cl_in.curl, col=cmb_col
;;          oplot, l_in, fl_in*abs(cl_in.G_T), col=cmb_col
;; ;         oplot, lx, fl_x*clxe, col=xcol[1]
;; ;         oplot, lx, fl_x*abs(clxte), col=xcol[3]
;; ;;         legendastro, ['CMB', $
;; ;;                       '!7D!3TT', $
;; ;;                       '!7D!3EE', $
;; ;;                       '!7D!3BB', $
;; ;;                       '!7D!3TE'], line=0, $
;; ;;                      col=[cmb_col, xcol[0:3]]
;;          plot_init = 1
;;       endif
;;                                 ;oplot, lx, fl_x*clxb, col=xcol[2]
;;       oplot, lx, fl_x*clxb, col=70 ; col=mycol[i_eta]
;;       w = where( lx eq 1000)
;;       xyouts, 1000, fl_x[w]*clxb[w]*1.1, "!7g!3 = 10!u"+strtrim( long(eta_pow[i_eta]),2)+"!n", col=70 ;, col=mycol[i_eta]
;; 
;;       wtol = where( l_in eq 90)
;;       ;; epsilon_list[i_eta] = sqrt( tol*cl_in[wtol].curl/clxb[w])
;;       ;; epsilon_list[i_eta] = sqrt( cosmic_var[wtol]/clxb[w])
;; 
;;       cosmic_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.temperature
;;       epsilon_list[i_eta,0] = sqrt( cosmic_var[wtol]/clxt[w])
;; 
;;       cosmic_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.gradient
;;       epsilon_list[i_eta,1] = sqrt( cosmic_var[wtol]/clxe[w])
;; 
;;       cosmic_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.curl
;;       epsilon_list[i_eta,2] = sqrt( cosmic_var[wtol]/clxb[w])
;; 
;;       cosmic_var = sqrt(1./(fsky*(2*l_in+1.d0)))*sqrt(cl_in.g_t^2 + cl_in.temperature*cl_in.gradient)
;;       epsilon_list[i_eta,3] = sqrt( cosmic_var[wtol]/abs(clxte[w]))
;; 
;;       ;; Output ascii file
;;       printf, lu, strtrim(nu_list[inu],2)+", 10^"+$
;;               strtrim(eta_pow[i_eta],2)+", "+$
;;               strtrim(epsilon_list[i_eta,0],2)+", "+$
;;               strtrim(epsilon_list[i_eta,1],2)+", "+$
;;               strtrim(epsilon_list[i_eta,2],2)+", "+$
;;               strtrim(epsilon_list[i_eta,3],2)
;; 
;;       legendastro, [strtrim(long(nu_list[inu]),2)+' GHz'], /right
;;       outplot, /close, /verb
;;    endfor
;; ;   outplot, /close, /verb
;; endfor
;; close, lu
;; free_lun, lu
;; 
;; ;; A single eta and display all freqs
;; outplot, file='fg_res_spectra_all_nu', png=png, ps=ps
;; eta_pow = -3
;; eta = 10.d0^eta_pow
;; make_ct, n_nu, my_col
;; for inu=0, n_nu-1 do begin
;; 
;;    ;; Scale foregrounds to the reference frequencies
;;    i_dust = i_dust_ref * (nu_list[inu]/545.d0)^beta_dust
;;    q_dust = q_dust_ref * (nu_list[inu]/353.d0)^beta_dust
;;    u_dust = u_dust_ref * (nu_list[inu]/353.d0)^beta_dust
;; 
;;    i_sync = i_sync_ref * (nu_list[inu]/0.408)^beta_sync
;;    q_sync = q_sync_ref * (nu_list[inu]/30.d0)^beta_sync
;;    u_sync = u_sync_ref * (nu_list[inu]/30.d0)^beta_sync
;; 
;;    xmap[*,0] = cmb_maps[*,0]^2 + $
;;                2.*cmb_maps[*,0]*eta*(i_dust + i_sync) + $
;;                eta^2*(i_dust + i_sync)^2
;;    xmap[*,1] = 2.d0*(cmb_maps[*,0]*cmb_maps[*,1] + $
;;                      eta*(q_dust+q_sync)*cmb_maps[*,0] + $
;;                      eta*(i_dust+i_sync)*cmb_maps[*,1])
;;    xmap[*,2] = 2.d0*(cmb_maps[*,0]*cmb_maps[*,2] + $
;;                      eta*(u_dust+u_sync)*cmb_maps[*,0] + $
;;                      eta*(i_dust+i_sync)*cmb_maps[*,2])
;;    
;;    ispice, xmap, cover, lx, clxt, clxe, clxb, clxte
;;    lx    = lx[2:*]
;;    clxt  = clxt[2:*]
;;    clxe  = clxe[2:*]
;;    clxb  = clxb[2:*]
;;    clxte = clxte[2:*]
;;    
;;    ;; Pixel window functions
;;    wlpix = mrdfits( !nika.soft_dir+"/idl-libs/Healpix/data/pixel_window_n"+$
;;                     string( nside, format='(I4.4)')+".fits", 1)
;;    ;; remove monopole and dipole
;;    wlpix = wlpix[2:ncl_x-1]
;;    
;;    ;; Derive spec
;;    w100 = where( l_in eq 100)
;;    wx100 = where( lx eq 100)
;;    
;;    if inu eq 0 then begin
;;       yra = [1d-7, 1d15]
;;       xra = [2, 3000]
;;       thick=2
;;       if ps eq 0 then wind, 1, 1, /free, /large
;;       plot, lx, fl_x*clxt, /xs, /ys, xra=xra, yra=yra, $
;;             xtitle='!12l', $
;;             ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n', $
;;             /ylog, /xlog, /nodata
;;       xyouts, 2.2, (fl_in*cl_in.temperature)[2]*1.3, 'TT'
;;       xyouts, 2.2, (fl_in*cl_in.gradient)[2]*1.3, 'EE'
;;       xyouts, 2.2, (fl_in*cl_in.curl)[2]*1.3, 'BB (r=0.001)'
;;       xyouts, 2.2, (fl_in*abs(cl_in.g_t))[2]*1.3, 'TE'
;;       oplot, l_in, fl_in*cl_in.temperature, col=0
;;       oplot, l_in, fl_in*cl_in.gradient, col=cmb_col
;;       oplot, l_in, fl_in*cl_in.curl, col=cmb_col
;;       oplot, l_in, fl_in*abs(cl_in.G_T), col=cmb_col
;;       legendastro, '!7g!3 = 10!u'+strtrim(eta_pow,2)+'!n', /right
;;       legendastro, ['CMB', $
;;                     '!7D!3BB ('+string(nu_list[inu],form='(I3.3)')+" GHz)"], $
;;                    line=0, col=[cmb_col, my_col]
;;    endif
;;    oplot, lx, fl_x*clxb, col=my_col[inu]
;; endfor
;; outplot, /close, /verb
;; stop
;; 
;; ;;-----------------------------------------------------
;; ;; Pure NL
;; yra = [1d-8, 1d5]
;; nl_res_file = 'nl_res_nl.dat'
;; eps_pow = [-2,-3,-4,-5,-6,-7,-8]
;; eps_list = 10.d0^eps_pow
;; n_eps = n_elements(eps_list)
;; openw, lu, nl_res_file, /get_lun
;; printf, lu, "#nu (GHz), epsilon, epsilonT, epsilonE, epsilonB, epsilonTE"
;; for inu=0, n_nu-1 do begin
;;    plot_init = 0
;; 
;;    ;; Scale foregrounds to the reference frequencies
;;    i_dust = i_dust_ref * (nu_list[inu]/545.d0)^beta_dust
;;    q_dust = q_dust_ref * (nu_list[inu]/353.d0)^beta_dust
;;    u_dust = u_dust_ref * (nu_list[inu]/353.d0)^beta_dust
;; 
;;    i_sync = i_sync_ref * (nu_list[inu]/0.408)^beta_sync
;;    q_sync = q_sync_ref * (nu_list[inu]/30.d0)^beta_sync
;;    u_sync = u_sync_ref * (nu_list[inu]/30.d0)^beta_sync
;; 
;;    ;; Loop on residuals amplitude
;;    for ieps=0, n_eps-1 do begin
;;       outplot, file='nl_res_spectra_nu'+strtrim( long(nu_list[inu]),2)+$
;;                'eta10'+strtrim(eps_pow[ieps],2), png=png, ps=ps
;; 
;;       i_total = cmb_maps[*,0] + i_dust + i_sync
;;       q_total = cmb_maps[*,1] + q_dust + q_sync
;;       u_total = cmb_maps[*,2] + u_dust + u_sync
;; 
;;       xmap[*,0] = cmb_maps[*,0]^2 + $
;;                   2.*cmb_maps[*,0]*eta*(i_dust + i_sync) + $
;;                   eta^2*(i_dust + i_sync)^2
;;       xmap[*,1] = 2.d0*(cmb_maps[*,0]*cmb_maps[*,1] + $
;;                         eta*(q_dust+q_sync)*cmb_maps[*,0] + $
;;                         eta*(i_dust+i_sync)*cmb_maps[*,1])
;;       xmap[*,2] = 2.d0*(cmb_maps[*,0]*cmb_maps[*,2] + $
;;                         eta*(u_dust+u_sync)*cmb_maps[*,0] + $
;;                         eta*(i_dust+i_sync)*cmb_maps[*,2])
;; 
;;       ispice, xmap, cover, lx, clxt, clxe, clxb, clxte
;;       lx    = lx[2:*]
;;       clxt  = clxt[2:*]
;;       clxe  = clxe[2:*]
;;       clxb  = clxb[2:*]
;;       clxte = clxte[2:*]
;;       
;;       ;; Pixel window functions
;;       wlpix = mrdfits( !nika.soft_dir+"/idl-libs/Healpix/data/pixel_window_n"+$
;;                        string( nside, format='(I4.4)')+".fits", 1)
;;       ;; remove monopole and dipole
;;       wlpix = wlpix[2:ncl_x-1]
;;       
;;       ;; Derive spec
;;       w100 = where( l_in eq 100)
;;       wx100 = where( lx eq 100)
;;       
;;       if plot_init eq 0 then begin
;;          ;; All spectra on the same plot
;;          xra = [2, 3000]
;;          thick=2
;;          if ps eq 0 then wind, 1, 1, /free, /large
;;          plot, lx, fl_x*clxt, /xs, /ys, xra=xra, yra=yra, $
;;                xtitle='!12l', $
;;                ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n', $
;;                /ylog, /xlog, /nodata
;;          xyouts, 2.2, (fl_in*cl_in.temperature)[2]*1.3, 'TT'
;;          xyouts, 2.2, (fl_in*cl_in.gradient)[2]*1.3, 'EE'
;;          xyouts, 2.2, (fl_in*cl_in.curl)[2]*1.3, 'BB (r=0.001)'
;;          xyouts, 2.2, (fl_in*abs(cl_in.g_t))[2]*1.3, 'TE'
;; ;         oplot, lx, fl_x*clxt, col=xcol[0]
;;          oplot, l_in, fl_in*cl_in.temperature, col=0
;;          oplot, l_in, fl_in*cl_in.gradient, col=cmb_col
;;          oplot, l_in, fl_in*cl_in.curl, col=cmb_col
;;          oplot, l_in, fl_in*abs(cl_in.G_T), col=cmb_col
;; ;         oplot, lx, fl_x*clxe, col=xcol[1]
;; ;         oplot, lx, fl_x*abs(clxte), col=xcol[3]
;;          legendastro, ['CMB', $
;;                        '!7D!3TT', $
;;                        '!7D!3EE', $
;;                        '!7D!3BB', $
;;                        '!7D!3TE'], line=0, $
;;                       col=[cmb_col, xcol[0:3]]
;;          xyouts, 50, 1, 'pure NL', orient=45, charsize=2
;;          plot_init = 1
;;       endif
;;                                 ;oplot, lx, fl_x*clxb, col=xcol[2]
;;       oplot, lx, fl_x*clxb, col=70 ; col=mycol[ieps]
;;       w = where( lx eq 1000)
;;       xyouts, 1000, fl_x[w]*clxb[w]*1.1, "!7e!3 = 10!u"+strtrim( long(eps_pow[ieps]),2)+"!n", col=70;, col=mycol[ieps]
;; 
;;       wtol = where( l_in eq 90)
;;       ;; epsilon_list[ieps] = sqrt( tol*cl_in[wtol].curl/clxb[w])
;;       ;; epsilon_list[ieps] = sqrt( cosmic_var[wtol]/clxb[w])
;; 
;;       cosmic_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.temperature
;;       epsilon_list[ieps,0] = sqrt( cosmic_var[wtol]/clxt[w])
;; 
;;       cosmic_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.gradient
;;       epsilon_list[ieps,1] = sqrt( cosmic_var[wtol]/clxe[w])
;; 
;;       cosmic_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.curl
;;       epsilon_list[ieps,2] = sqrt( cosmic_var[wtol]/clxb[w])
;; 
;;       cosmic_var = sqrt(1./(fsky*(2*l_in+1.d0)))*sqrt(cl_in.g_t^2 + cl_in.temperature*cl_in.gradient)
;;       epsilon_list[ieps,3] = sqrt( cosmic_var[wtol]/abs(clxte[w]))
;; 
;;       ;; Output ascii file
;;       printf, lu, strtrim(nu_list[inu],2)+", 10^"+$
;;               strtrim(eps_pow[ieps],2)+", "+$
;;               strtrim(epsilon_list[ieps,0],2)+", "+$
;;               strtrim(epsilon_list[ieps,1],2)+", "+$
;;               strtrim(epsilon_list[ieps,2],2)+", "+$
;;               strtrim(epsilon_list[ieps,3],2)
;; 
;;       outplot, /close, /verb
;;    endfor
;; endfor
;; close, lu
;; free_lun, lu
;; 
;; 
;; 
;; 
;; 
;; 
;; 
;; 
;; 
;; 
;; spawn, "cat "+fg_res_file
;; 
;; readcol, fg_res_file, nu, eta_test, epst, epse, epsb, epste, $
;;          format='A,A,A,A,A,A', delim=',', comment='#'
;; all_eps = dblarr(4,n_elements(nu))
;; all_eps[0,*] = double(epst)
;; all_eps[1,*] = double(epse)
;; all_eps[2,*] = double(epsb)
;; all_eps[3,*] = double(epste)
;; spec_title = ['T', 'E', 'B', 'TE']
;; make_ct, n_nu, ct
;; if ps eq 0 then wind, 1, 1, /free, /large
;; my_multiplot, 2, 2, pp, pp1, /rev
;; outplot, file='fg_res_nl_plot', ps=ps, png=png
;; xra = minmax(eta_list)*[0.1,10]
;; yra = minmax(all_eps)
;; for ispec=0, 3 do begin
;;    for i_nu=0, n_nu-1 do begin
;;       w = where( double(nu) eq nu_list[i_nu], nw)
;;       if i_nu eq 0 then begin
;;          plot, eta_list, all_eps[ispec,w], /xlog, /ylog, xra=xra, yra=yra, $
;;                xtitle='!7g!3', ytitle='!7e!3 ('+strtrim(spec_title[ispec],2)+')', $
;;                position=pp1[ispec,*], /noerase
;;          legendastro, strtrim( long(nu_list),2)+" GHz", col=ct, line=0, /right
;;       endif
;;       oplot, eta_list, all_eps[ispec,w], psym=-8, col=ct[i_nu]
;;    endfor
;; endfor
;; outplot, /close, /verb

end
