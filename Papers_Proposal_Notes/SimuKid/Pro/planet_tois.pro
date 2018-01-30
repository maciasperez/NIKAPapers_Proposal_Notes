

;; Conversion Jy to Hz based on N2R11
jy2hz = 1500./30

col_rf = 70
col_cf = 250

png = 0 ; 1
ps  = 1
nside = 2048

;; We first fix the HWP rotation frequency
fwhm = 11. ; NIKA2
hwp_rot_freq = 3.d0 ; 3.d0 like NIKA
n_harmonics = 8
f_sampling = long( (12.*hwp_rot_freq) > (2*n_harmonics*hwp_rot_freq))

npts_per_fwhm_list = [3, 5]
nspeed = n_elements( npts_per_fwhm_list)
my_multiplot, 2, 2, pp, pp1, /rev
if ps eq 0 then wind, 1, 1, /free, /large

;; premiere fois qu'on definit kid_model, donc on l'initialise,
;; structure kid_model
kid_model = {z0:1.d0, x1:1.d0, x2:1.d0, Qe:1.d0, Qi:1.d0, f0:1.d0}
kid_model.x1 = 5d-1              ;3d0 ; 0.d0
kid_model.x2 = 5d-1              ;3d0 ; 0.d0
kid_model.z0 = 50.d0
kid_model.Qe = 1d4              ;2d4 ;5.2d4
kid_model.Qi = 2d4              ;5.2d4 ; 1.9d5
kid_model.f0 = 1.8d9            ;1.273d9 ; Hz
delta_f      = 5d3              ;2d3 ; Martino, or 4, 5, 8

;; Frequences de lecture avec la modulation
f_m = kid_model.f0 - delta_f/2.
f_p = kid_model.f0 + delta_f/2.

;; Oversample to prepare I, Q, dI, dQ
fsampl_factor = 20L
f_hf = fsampl_factor*long(f_sampling) ;f_sampling*40.d0

;; Time parameters HF and LF
obs_time_hour = 1./60d0
nsn = round(obs_time_hour*3600.d0*f_hf)
time_hf = dindgen(nsn)/f_hf

gplanet_speed = dblarr(nspeed,nsn)
for ispeed=0, nspeed-1 do begin
   npts_per_fwhm_min = npts_per_fwhm_list[ispeed]
   vmax = fwhm/npts_per_fwhm_min * 4 * hwp_rot_freq ; arcsec/s

;; To produce signal timelines from pixelized maps 
f_nyquist = 1./(fwhm/3./vmax)
nsn_nyquist = round( obs_time_hour*3600.d0*f_nyquist)

;; Planet timeline
   scan_speed = vmax            ; arcmin/s
   sigma = fwhm*!fwhm2sigma
   
   x_hf = scan_speed*time_hf
   x0 = x_hf[nsn/4]
   x_hf_min = min(x_hf)
   t_planet_hf = time_hf[nsn/4]
   
;; Planet gaussian profile
   gplanet = exp(-(x_hf-x0)^2/(2.*sigma^2))
   gplanet_speed[ispeed,*] = gplanet
   
;; Main loop
   in_flux_list = [1, 100, 1000, 2000] ; Jy
   nflux = n_elements(in_flux_list)
   for iflux=0, nflux-1 do begin
      toi_planet_jy = gplanet * in_flux_list[iflux]
      freso_planet  = -toi_planet_jy*jy2hz
      
      freso2toi, freso_planet, kid_model, delta_f, toi_rf_planet_hz, toi_cf_planet_hz, $
                 i=i, q=q, stop=stop, xc=xc, yc=yc, r=r
      if defined(toi_rf_res_jy) eq 0 then begin
         ;; Define the ouput LF time vector:
         ;; need to add 0.5/f_sampling to compute the time corresponding to the average
         ;; of the 20-40 points interval taken in iqdidq to produce a final sample from
         ;; the "oversampled" i,q timelines.
         ni = n_elements(toi_rf_planet_hz)
         time_lf = dindgen(ni)/f_sampling + 0.5/f_sampling
         
         toi_rf_res_jy = dblarr(nspeed,nflux,ni)
         toi_cf_res_jy = dblarr(nspeed,nflux,ni)
      endif
      
      ;; convert back to Jy for the plot
      toi_rf_planet_jy = toi_rf_planet_hz/jy2hz
      toi_cf_planet_jy = toi_cf_planet_hz/jy2hz
      
      toi_rf_res_jy[ispeed,iflux,*] = toi_rf_planet_jy
      toi_cf_res_jy[ispeed,iflux,*] = toi_cf_planet_jy
   endfor
endfor


;;-------------------------------------
;; Plot of timelines
;make_ct, nflux, flux_col
flux_col = [70, 150, 200, 250]

xra = [14,16]
symlist = [1,4,5,6]
outplot, file='planet_profiles', png=png, ps=ps
for ispeed=0, nspeed-1 do begin
   plot, time_hf, gplanet_speed[ispeed,*], xra=xra, yra=yra, /xs, /ys, $
         position=pp[0,ispeed,*], xtitle='time (sec)', ytitle='Output/input Flux', /noerase
   for iflux=0, nflux-1 do begin
      norm = in_flux_list[iflux]
      oplot, time_lf, toi_rf_res_jy[ispeed,iflux,*]/norm, $
             col=flux_col[iflux], psym=symlist[iflux]
   endfor
   legendastro, strtrim(in_flux_list,2)+" Jy", col=flux_col, psym=symlist
   legendastro, ['RfdIdQ', strtrim(npts_per_fwhm_list[ispeed],2)+" pts/FWHM"], /right
   plot, time_hf, gplanet_speed[ispeed,*], xra=xra, yra=yra, /xs, /ys, $
         position=pp[1,ispeed,*], /noerase, xtitle='time (sec)', $
         ycharsize=1d-10
   for iflux=0, nflux-1 do begin
      norm = in_flux_list[iflux]
      oplot, time_lf, toi_cf_res_jy[ispeed,iflux,*]/norm, $
             col=flux_col[iflux], psym=symlist[iflux]
   endfor
   legendastro, strtrim(in_flux_list,2)+" Jy", col=flux_col, psym=symlist
   legendastro, ['Cf', strtrim(npts_per_fwhm_list[ispeed],2)+" pts/FWHM"], /right
endfor
outplot, /close, /verb

;; 
;; 
;;    xra = t_planet_hf + [-1, 1]*0.5
;;    norm_ampl = max(toi_planet_jy)
;;    if did_plot eq 0 then begin
;;       plot,  time_hf, toi_planet_jy/norm_ampl, /xs, position=pp[0,0,*], /noerase, $
;;              xra=xra, xtitle='Time (sec)', ytitle='Output Flux / Input flux'
;;       did_plot = 1
;;    endif
;; 
;;    wfit = where( abs(time_lf-t_planet_hf) le 2,nwfit)
;;    junk = gaussfit( time_lf[wfit], toi_rf_planet_jy[wfit], a, nterms=5)
;;    flux = a[0]
;;    oplot, time_lf, toi_rf_planet_jy/a[0], psym=8, syms=0.5, col=flux_col[iflux]
;; 
;;    legendastro, strtrim(in_flux_list,2)+" Jy", line=0, col=flux_col
;; endfor
;; xyouts, 14.8, 0.4, orien=45, 'WHY are the output FWHM different ?!'

end
