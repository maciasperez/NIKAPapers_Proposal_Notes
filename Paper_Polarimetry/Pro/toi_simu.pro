
;; Quick simulation of polarized point source observation to explain the map
;; making principle

nk_default_param, scan_params
nk_default_info, info

;; Scanning strategy
az_width = 6.5*60+2*60. ; large to have all kid off the source, but completely off source some fraction of the time

scan_params.decor_cm_dmin = 10.

az_width = 6.5*60-10
az_speed = 25.
el_width = 50
el_step = 2.
n_subscans = round(el_width/el_step)
scan_params = create_struct( scan_params, $
                             "f_sampling" , 47.68d0, $ ; Hz
                             "az_width"   , az_width, $ ; 400.d0                 ; arcsec
                             "el_step"    , el_step, $   ; arcsec
                             "angle_deg"  , 0.d0, $ ; degrees
                             "az_speed"   , az_speed, $   ; arcsec/s
                             "el_speed"   , 60.d0, $   ; arcsec/s, intersubscan, place holder
                             "n_subscans" , n_subscans)

scan_params.map_xsize = (3*az_width) > 700
scan_params.map_ysize = scan_params.map_xsize
scan_params.map_reso = 2.d0
scan_params.map_proj = 'radec'
nk_init_grid, scan_params, info, input_grid

nu_hwp = 2.98d0


;;----------------------------------------------------------

;; Scanning strategy
el_min = -scan_params.n_subscans/2*scan_params.el_step
n_samples_per_subscan     = round(scan_params.az_width/scan_params.az_speed*scan_params.f_sampling)
n_sample_per_intersubscan = round(scan_params.el_step/scan_params.el_speed*scan_params.f_sampling)
nsn = scan_params.n_subscans*n_samples_per_subscan
ofs_az = dblarr(nsn)
ofs_el = dblarr(nsn)
scan_obs_time = nsn/scan_params.f_sampling

;; for i=0, scan_params.n_subscans-1 do begin
;;    ofs_az[i*n_samples_per_subscan: (i+1)*n_samples_per_subscan-1] = $
;;       -scan_params.az_width/2.d0 + dindgen( n_samples_per_subscan)/scan_params.f_sampling*scan_params.az_speed
;;    ofs_el[i*n_samples_per_subscan: (i+1)*n_samples_per_subscan-1] = $
;;       el_min + i*scan_params.el_step
;; endfor

last_az = -scan_params.az_width/2.d0
ofs_az = [0]
ofs_el = [0]
for i=0, scan_params.n_subscans-1 do begin
   sign = 1
   if (i mod 2) eq 1 then sign=-1
   az = sign*dindgen( n_samples_per_subscan)/scan_params.f_sampling*scan_params.az_speed
   ofs_az = [ofs_az, az + last_az]
   last_az = last_az + az[n_samples_per_subscan-1]
   ofs_el = [ofs_el, dblarr(n_samples_per_subscan) + el_min + i*scan_params.el_step]
   last_el = el_min + i*scan_params.el_step

   ;; intersubscan
   el = dindgen(4)/3*scan_params.el_step
   ofs_az = [ofs_az, dblarr(4)+last_az]
   ofs_el = [ofs_el, last_el + el]
endfor

ofs_az = ofs_az[1:*]
ofs_el = ofs_el[1:*]

;; Compute pointing per kid and project
kidpar = mrdfits( !nika.off_proc_dir+"/kidpar_20160303s89_90_91_noskydip_v2.fits", 1, /silent)
nkids = n_elements(kidpar)

;; get_kid_pointing
alpha = dblarr(nsn) + 45*!dtor ; alpha_nasmyth( data.el)
dx  = kidpar.nas_x - kidpar.nas_center_x
dy  = kidpar.nas_y - kidpar.nas_center_y
daz = cos(alpha)##dx - sin(alpha)##dy
del = sin(alpha)##dx + cos(alpha)##dy

dra  = -daz + ofs_az##( dblarr(nkids)+1)
ddec = -del + ofs_el##( dblarr(nkids)+1)

time = dindgen(nsn)/scan_params.f_sampling
omega = 2d0*!dpi*nu_hwp*time mod (2.d0*!dpi)

;; Signal
Iflux = 10
Qflux = 1
Uflux = 0
fwhm = 11.d0

sigma_beam = fwhm*!fwhm2sigma
ikid = where( kidpar.numdet eq 860)
s_toi_i = Iflux * exp(-(dra[ikid,*]^2+ddec[ikid,*]^2)/(2.*sigma_beam^2))
s_toi_p = Qflux * exp(-(dra[ikid,*]^2+ddec[ikid,*]^2)/(2.*sigma_beam^2))*cos(4.d0*omega) + $
          Uflux * exp(-(dra[ikid,*]^2+ddec[ikid,*]^2)/(2.*sigma_beam^2))*sin(4.d0*omega)

s_toi_i = reform(s_toi_i)
s_toi_p = reform(s_toi_p)


sigma_wn = 0.01
noise = randomn( seed, nsn)*sigma_wn
toi = s_toi_i + s_toi_p + noise

power_spec,   toi,   scan_params.f_sampling, pwi, freq
power_spec, s_toi_i, scan_params.f_sampling, pwi_s
power_spec, s_toi_p, scan_params.f_sampling, pwi_p

power_spec,     toi*cos(4.d0*omega), scan_params.f_sampling, pwq
power_spec, s_toi_p*cos(4.d0*omega), scan_params.f_sampling, pwq_p
power_spec, s_toi_i*cos(4.d0*omega), scan_params.f_sampling, pwq_i


nf    = 1000
f_max = 30                      ; Hz
f     = dindgen( nf)/(nf-1)*f_max
f[0]  = 1e-10 ; for a nice disply on log plots
sigma_t = fwhm*!fwhm2sigma/az_speed
sigma_k = 1.0d0/(2.0d0*!dpi*sigma_t)

;; because of the sim sampling, the effective max of the TOI is not exactly the
;; nominal fluxes
Ia = max(pwi)
Qa = max(pwi[where( abs(freq-4*nu_hwp) lt 0.1)])


ps =1 
;; outplot, file='../figures/toi_simu_pol', png=png, ps=ps
;; plot_oo, freq, pwi, /xs, ys=5, $
;;          yra=yra, xra=xra, position=pp1[0,*], xchars=1e-10, chars=1.2, charthick=3
;; oplot, freq, pwi_s, col=100
;; oplot, freq, pwi_p, col=150
;; outplot, /close
;; pstopdf_crop, '../figures/toi_simu_pol'
;; stop
;; plot en log-lin
;; if defined(ps) eq 0 then wind, 1, 1, /free, /large
;; outplot, file='../figures/toi_simu', png=png, ps=ps
set_plot, 'ps'
device, filename='../figures/toi_simu.ps',  bits_per_pixel=256,/color
  
my_multiplot, 2, 2, pp, pp1, /full, /dry, ymargin=0.1, /rev
x1 = 6
xra = [1e-3, x1]
yra = [1e-8, 1e2]
plot_oo, freq, pwi, /xs, ys=5, $
         yra=yra, xra=xra, position=pp1[0,*], xchars=1e-10, chars=1.2, charthick=3
oplot, freq, pwi_s, col=100
oplot, freq, pwi_p, col=150
oplot, f, Ia*exp(-f^2/(2.0d0*sigma_k^2)), col=250
oplot, f, Qa*exp(-(4.d0*nu_hwp-f)^2/(2.0d0*sigma_k^2)), col=250
for i=1, 2 do begin
   oplot, [1,1]*i*nu_hwp, [1e-20, 1e10], line=2
   xyouts, i*nu_hwp*1.05, 1e1, strtrim(i,2)+"!7m!3!dP!n", /data, chars=1.2, charthick=3
endfor
legendastro, ['Signal+Noise', 'Signal I', 'Signal P', 'Beam cut off'], $
            line=0, col=[!p.color, 100, 150, 250], box=0, /bottom, chars=1.2, charthick=3
legendastro, 'Total power', box=0, chars=1.2, charthick=3
axis, yaxis=0, yrange=yra, /ylog, /ys, ytitle='Arb. Units x Hz!u-1/2!n', chars=1.2, charthick=3

plot, freq, pwi, xra=[x1, max(freq)], /xs, ys=5, yra=yra, /ylog, $
      position=pp1[1,*], /noerase, xchars=1e-10, chars=1.2, charthick=3
oplot, freq, pwi_s, col=100
oplot, freq, pwi_p, col=150
oplot, f, Ia*exp(-f^2/(2.0d0*sigma_k^2)), col=250
oplot, f, Qa*exp(-(4.d0*nu_hwp-f)^2/(2.0d0*sigma_k^2)), col=250
for i=3, 7 do begin
   oplot, [1,1]*i*nu_hwp, [1e-20, 1e10], col=0, line=2
   xyouts, i*nu_hwp+.1, 1e1, strtrim(i,2)+"!7m!3!dP!n", col=0, /data,  chars=1.2, charthick=3
endfor
axis, yaxis=1, yrange=yra, /ylog, /ys, ycharsize=1e-10, chars=1.2, charthick=3

plot_oo, freq, pwq, /xs, ys=5, $
         yra=yra, xra=xra, position=pp1[2,*], /noerase, chars=1.2, charthick=3
oplot, freq, pwq_i, col=100
oplot, freq, pwq_p, col=150
oplot, f, Qa*exp(-f^2/(2.0d0*sigma_k^2)), col=250
for i=1, 2 do begin
   oplot, [1,1]*i*nu_hwp, [1e-20, 1e10], col=0, line=2
endfor
legendastro, 'Demodulated', box=0, chars=1.2, charthick=3
axis, yaxis=0, yrange=yra, /ylog, /ys, ytitle='Arb. Units x Hz!u-1/2!n', chars=1.2, charthick=3
xyouts, 0.44, 0.02, 'Frequency [Hz]', /norm,  chars=1.2, charthick=3

plot, freq, pwq, xra=[x1, max(freq)], /xs, ys=5, yra=yra, /ylog, $
      position=pp1[3,*], /noerase,  chars=1.2, charthick=3
oplot, freq, pwq_p, col=150
oplot, freq, pwq_i, col=100
oplot, f, Qa*exp(-f^2/(2.0d0*sigma_k^2)), col=250
for i=3, 7 do begin
   oplot, [1,1]*i*nu_hwp, [1e-20, 1e10], col=0, line=2
endfor
axis, yaxis=1, yrange=yra, /ylog, /ys, ycharsize=1e-10, chars=1.2, charthick=3
;; outplot, /close
device, /close
set_plot, 'x'
pstopdf_crop, '../figures/toi_simu'
my_multiplot, /res

;; outplot, file='../figures/toi_simu', png=png, ps=ps
;; my_multiplot, 2, 1, pp, pp1, /full, /dry, ymargin=0.1
;; x1 = 6
;; xra = [1e-3, x1]
;; yra = [1e-12, 1e2]
;; plot_oo, freq, pw_toi, /xs, ys=5, $
;;          yra=yra, xra=xra, position=pp1[0,*], chars=1.2, charthick=3
;; oplot, freq, pw_s, col=150
;; oplot, f, Ia*exp(-f^2/(2.0d0*sigma_k^2)), col=250
;; oplot, f, Qa*exp(-(4.d0*nu_hwp-f)^2/(2.0d0*sigma_k^2)), col=250
;; for i=1, 2 do begin
;;    oplot, [1,1]*i*nu_hwp, [1e-20, 1e10], col=70, line=2
;;    xyouts, i*nu_hwp*1.1, 1e1, strtrim(i,2)+"!7x!3", col=70, /data
;; endfor
;; legendastro, ['Signal+Noise', 'Signal', 'Beam cut off'], $
;;             line=0, col=[!p.color, 150, 250], box=0, /bottom
;; axis, yaxis=0, yrange=yra, /ylog, /ys, ytitle='Arb. Units x Hz!u-1/2!n', chars=1.2, charthick=3
;; xyouts, 0.5, 0.02, 'Frequency (Hz)', /norm, chars=1.2, charthick=3

;; plot, freq, pw_toi, xra=[x1, max(freq)], /xs, ys=5, yra=yra, /ylog, $
;;       position=pp1[1,*], /noerase, chars=1.2, charthick=3
;; oplot, freq, pw_s, col=150
;; oplot, f, Ia*exp(-f^2/(2.0d0*sigma_k^2)), col=250
;; oplot, f, Qa*exp(-(4.d0*nu_hwp-f)^2/(2.0d0*sigma_k^2)), col=250
;; for i=3, 7 do begin
;;    oplot, [1,1]*i*nu_hwp, [1e-20, 1e10], col=70, line=2
;;    xyouts, i*nu_hwp+.1, 1e1, strtrim(i,2)+"!7x!3", col=70, /data
;; endfor
;; axis, yaxis=1, yrange=yra, /ylog, /ys, ycharsize=1e-10
;; outplot, /close
;; my_multiplot, /res


end
