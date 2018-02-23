
;; From hwp_flux.pro
;;--------------------

;; Realizations
nmc = 30 ; 100

;; HWP power
hwp_max_ampl_jy = 1000 ; 100. ; 50.

;; Conversion Jy to Hz based on N2R11
jy2hz = 1500./30

col_rf = 70
col_cf = 250
png = 0 ; 1
ps  = 0

;; We first fix the HWP rotation frequency
fwhm = 11. ; NIKA2
hwp_rot_freq = 3.d0 ; 3.d0 like NIKA
n_harmonics = 8
f_sampling = long( (12.*hwp_rot_freq) > (2*n_harmonics*hwp_rot_freq))

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

;; Planet fluxes
nflux = 30
flux_min = 1.d0
flux_max = 500.d0

;;--------------------------------------------------------------
;; Scanning speeds
npts_per_fwhm_list = [3, 5]
nspeed = n_elements( npts_per_fwhm_list)

;; Planet fluxes
flux_list = dindgen(nflux)/(nflux-1.)*(flux_max-flux_min) + flux_min

;; Init result arrays
planet_rf_flux_results     = dblarr(nspeed, nmc, nflux)
planet_cf_flux_results     = dblarr(nspeed, nmc, nflux)
planet_hwp_rf_flux_results = dblarr(nspeed, nmc, nflux)
planet_hwp_cf_flux_results = dblarr(nspeed, nmc, nflux)

calib_planet_rf       = dblarr(nspeed,nmc)
epsilon_planet_rf     = dblarr(nspeed,nmc)
const_planet_rf       = dblarr(nspeed,nmc)
calib_planet_hwp_rf   = dblarr(nspeed,nmc)
epsilon_planet_hwp_rf = dblarr(nspeed,nmc)
const_planet_hwp_rf   = dblarr(nspeed,nmc)

calib_planet_cf       = dblarr(nspeed,nmc)
epsilon_planet_cf     = dblarr(nspeed,nmc)
const_planet_cf       = dblarr(nspeed,nmc)
calib_planet_hwp_cf   = dblarr(nspeed,nmc)
epsilon_planet_hwp_cf = dblarr(nspeed,nmc)
const_planet_hwp_cf   = dblarr(nspeed,nmc)

;; All random numbers at once to minimize correlations
xn = reform( randomn( seed, n_harmonics*2*nmc), 2, n_harmonics*nmc) * hwp_max_ampl_jy

;; Gaussian fit parameters
nterms = 3

;; Main loop: produced TOI's and fit output fluxes
for ispeed=0, nspeed-1 do begin
   npts_per_fwhm_min = npts_per_fwhm_list[ispeed]
   vmax = fwhm/npts_per_fwhm_min * 4 * hwp_rot_freq    ; arcsec/s
   scan_speed = vmax                                   ; arcmin/s
   sigma = fwhm*!fwhm2sigma
   
   x_hf = scan_speed*time_hf
   x0 = x_hf[nsn/4]
   x_hf_min = min(x_hf)
   t_planet_hf = time_hf[nsn/4]
   gplanet = exp(-(x_hf-x0)^2/(2.*sigma^2))
   
   ;; Loop on hwp template
   for imc=0, nmc-1 do begin
      percent_status, imc, nmc, 10

      ;; HWP template
      an = xn[ 0, imc*n_harmonics:(imc+1)*n_harmonics-1]
      bn = xn[ 1, imc*n_harmonics:(imc+1)*n_harmonics-1]
      hwp_beta_jy = time_hf*0.d0
      for n=0, n_harmonics-1 do begin
         hwp_beta_jy += an[n]*cos(n*2.d0*!dpi*hwp_rot_freq*time_hf)
         hwp_beta_jy += bn[n]*sin(n*2.d0*!dpi*hwp_rot_freq*time_hf)
      endfor
      
      for iflux=0, nflux-1 do begin

         ;; Planet only
         freso = -gplanet*flux_list[iflux]*jy2hz
         freso2toi, freso, kid_model, delta_f, toi_rf_hz, toi_cf_hz, $
                    npts_avg=npts_avg, i=i_hwp, q=q_hwp
         toi_rf_jy = toi_rf_hz/jy2hz
         toi_cf_jy = toi_cf_hz/jy2hz

         if defined(time_lf) eq 0 then begin
            ni = n_elements(toi_rf_jy)
            time_lf = dindgen(ni)/f_sampling + 0.5/f_sampling
            ;; for the gaussian fit of the flux
            wfit = where( abs(time_lf-t_planet_hf) le 2.d0, nwfit)
         endif

         junk = gaussfit( time_lf, toi_rf_jy, a, nterms=nterms)
         planet_rf_flux_results[ispeed,imc,iflux] = a[0]
         junk = gaussfit( time_lf, toi_cf_jy, a, nterms=nterms)
         planet_cf_flux_results[ispeed,imc,iflux] = a[0]

         ;; HWP only (for subtraction)
         freso_hwp  =  -hwp_beta_jy*jy2hz
         freso2toi, freso_hwp, kid_model, delta_f, toi_rf_hwp_hz, toi_cf_hwp_hz, $
                    i=i_hwp, q=q_hwp
         toi_rf_hwp_jy = toi_rf_hwp_hz/jy2hz
         toi_cf_hwp_jy = toi_cf_hwp_hz/jy2hz

         ;; Planet + HWP
         freso = -(gplanet*flux_list[iflux]+hwp_beta_jy) * jy2hz
         freso2toi, freso, kid_model, delta_f, toi_rf_hz, toi_cf_hz, $
                    npts_avg=npts_avg, i=i_hwp, q=q_hwp, di=di, dq=dq
         toi_rf_jy = toi_rf_hz/jy2hz
         toi_cf_jy = toi_cf_hz/jy2hz
         ;; Subtract the HWP timeline not to bias the fit
         junk = gaussfit( time_lf, (toi_rf_jy-toi_rf_hwp_jy), a, nterms=nterms)
         planet_hwp_rf_flux_results[ispeed,imc,iflux] = a[0] ; flux
         junk = gaussfit( time_lf, (toi_cf_jy-toi_cf_hwp_jy), a, nterms=nterms)
         planet_hwp_cf_flux_results[ispeed,imc,iflux] = a[0] ; flux
      endfor
   endfor
endfor

;; Fit linearity coeff's
;wind, 1, 1, /free, /large
my_multiplot, 1, nspeed, pp, pp1, /rev
for ispeed=0, nspeed-1 do begin
   for imc=0, nmc-1 do begin
      ;; Planet only
      fit = poly_fit( flux_list, planet_rf_flux_results[ispeed,imc,*], 2, /double, /yfit, status=status)
      calib_planet_rf[  ispeed,imc] = 1.d0/fit[1]
      epsilon_planet_rf[ispeed,imc] = fit[2]
      const_planet_rf[  ispeed,imc] = fit[0]
      fit = poly_fit( flux_list, planet_cf_flux_results[ispeed,imc,*], 2, /double, /yfit, status=status)
      calib_planet_cf[  ispeed,imc] = 1.d0/fit[1]
      epsilon_planet_cf[ispeed,imc] = fit[2]
      const_planet_cf[  ispeed,imc] = fit[0]

      ;; Planet and hwp
      fit = poly_fit( flux_list, planet_hwp_rf_flux_results[ispeed,imc,*], 2, /double, /yfit, status=status)
      calib_planet_hwp_rf[  ispeed,imc] = 1.d0/fit[1]
      epsilon_planet_hwp_rf[ispeed,imc] = fit[2]
      const_planet_hwp_rf[  ispeed,imc] = fit[0]
      fit = poly_fit( flux_list, planet_hwp_cf_flux_results[ispeed,imc,*], 2, /double, /yfit, status=status)
      calib_planet_hwp_cf[  ispeed,imc] = 1.d0/fit[1]
      epsilon_planet_hwp_cf[ispeed,imc] = fit[2]
      const_planet_hwp_cf[  ispeed,imc] = fit[0]

      ;; plot, flux_list, flux_list, /xs, /ys, $
      ;;       xtitle='Input flux', ytitle='Output flux', $
      ;;       position=pp1[ispeed,*], /noerase
      ;; oplot, flux_list, planet_rf_flux_results[ispeed,imc,*,1], psym=1, col=col_rf
      ;; oplot, flux_list, planet_cf_flux_results[ispeed,imc,*,1], psym=4, col=col_cf
      ;; legendastro, strtrim(npts_per_fwhm_list[ispeed],2)+" pts/FWHM"
      ;; 
      ;; print, ""
      ;; print, "imc: "+Strtrim(imc,2)
      ;; print, "epsilon_planet_rf:     "+strtrim(epsilon_planet_rf[    ispeed,imc],2)
      ;; print, "epsilon_planet_cf:     "+strtrim(epsilon_planet_cf[    ispeed,imc],2)
      ;; print, "epsilon_planet_hwp_rf: "+strtrim(epsilon_planet_hwp_rf[ispeed,imc],2)
      ;; print, "epsilon_planet_hwp_cf: "+strtrim(epsilon_planet_hwp_cf[ispeed,imc],2)

   endfor
endfor

;; Output vs input flux
yra = [0.9, 1.01]
xra = minmax(flux_list)
imc = 0
if ps eq 0 then wind, 1, 1, /free, /large
outplot, file='flux_out_vs_in', png=png, ps=ps
plot, flux_list, planet_rf_flux_results[0,imc,*,1]/flux_list, $
      /xs, /ys, xra=xra, yra=yra, xtitle='Input flux', $
      ytitle='Output flux / input flux', /nodata
oplot, xra, xra*0. + 1.d0
oplot, flux_list, planet_rf_flux_results[0,imc,*]/flux_list, col=col_rf, line=2
oplot, flux_list, planet_rf_flux_results[1,imc,*]/flux_list, col=col_rf
oplot, flux_list, planet_cf_flux_results[0,imc,*]/flux_list, col=col_cf, line=2
oplot, flux_list, planet_cf_flux_results[1,imc,*]/flux_list, col=col_cf
legendastro, ['3pts/FWHM', '5pts/FWHM'], line=[2,0], /right
outplot, /close, /verb
stop

;; Distribution of non linearity coeffs
if ps eq 0 then wind, 1, 1, /free, /large
outplot, file='histos_epsilon', png=png, ps=ps
my_multiplot, 1, 1, ntot=nspeed*2, pp, pp1, /rev
for ispeed=0, nspeed-1 do begin
   np_histo, reform( epsilon_planet_hwp_rf[ispeed,*], nmc), $
             position=pp[ispeed,0,*], title='Epsilon planet_hwp RF', $
             /fit, /noerase, /fill, fcol=col_Rf
   legendastro, ['Max. HWP Ampl: '+strtrim(hwp_max_ampl_jy,2)+" Jy", $
                 strtrim(npts_per_fwhm_list[ispeed],2)+" pts/FWHM"]
   
   np_histo, reform( epsilon_planet_hwp_cf[ispeed,*],nmc), $
             position=pp[ispeed,1,*], title='Epsilon planet_hwp CF', $
             /fit, /noerase, /fill, fcol=col_cf
   legendastro, ['Max. HWP Ampl: '+strtrim(hwp_max_ampl_jy,2)+" Jy", $
                 strtrim(npts_per_fwhm_list[ispeed],2)+" pts/FWHM"]
endfor
outplot, /close, /verb


np_histo, reform( calib_planet_hwp_rf[ispeed,*], nmc), $
          position=pp1[0,*], title='Calib planet_hwp RF', /fit, /fill, fcol=col_rf
legendastro, ['Max. HWP Ampl: '+strtrim(hwp_max_ampl_jy,2)+" Jy"]

np_histo, reform( calib_planet_hwp_cf[ispeed,*], nmc), $
          position=pp1[1,*], title='Calib planet_hwp CF', /fit, /noerase, /fill, fcol=col_cf
legendastro, ['Max. HWP Ampl: '+strtrim(hwp_max_ampl_jy,2)+" Jy"]



end
