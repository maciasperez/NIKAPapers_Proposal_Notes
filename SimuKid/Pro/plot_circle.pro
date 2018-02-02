

;; Conversion Jy to Hz based on N2R11
jy2hz = 1500./30

col_rf = 70
col_cf = 250

png = 0 ; 1
ps  = 0
nside = 2048

npts_per_fwhm_min = 5 ; 3

;; We first fix the HWP rotation frequency
fwhm = 11. ; NIKA2
hwp_rot_freq = 3.d0 ; 3.d0 like NIKA
n_harmonics = 8
f_sampling = long( (12.*hwp_rot_freq) > (2*n_harmonics*hwp_rot_freq))
vmax = fwhm/npts_per_fwhm_min * 4 * hwp_rot_freq ; arcsec/s

;; To produce signal timelines from pixelized maps 
f_nyquist = 1./(fwhm/3./vmax)

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
nsn_nyquist = round( obs_time_hour*3600.d0*f_nyquist)

;; Planet timeline
scan_speed = vmax               ; arcmin/s
sigma = fwhm*!fwhm2sigma

x_hf = scan_speed*time_hf
x0 = x_hf[nsn/4]
x_hf_min = min(x_hf)
t_planet_hf = time_hf[nsn/4]

;; Planet gaussian profile
gplanet = exp(-(x_hf-x0)^2/(2.*sigma^2))

;; Main loop
;; in_flux_list = [1, 1000] ; Jy
in_flux_list = [100, 300, 700, 1000];, 10000, 1e10] ; Jy
nflux = n_elements(in_flux_list)

;wind, 1, 1, /free, /large
;my_multiplot, 2, 1, pp, pp1, /rev

; loop in reverse order to init plot ranges
delvarx, ires, qres
flux_col = [70, 100, 150, 250] ; dindgen(nflux)/(nflux-1)*(250-40) + 40
did_plot = 0
for iflux=nflux-1, 0, -1 do begin
   toi_planet_jy = gplanet * in_flux_list[iflux]
   freso_planet  = -toi_planet_jy*jy2hz
   freso2toi, freso_planet, kid_model, delta_f, toi_rf_planet_hz, toi_cf_planet_hz, $
              i=i, q=q, stop=stop, xc=xc, yc=yc, r=r
   
   ;; Rotate it for the picture
   alpha = 30*!dtor
   i1 = cos(alpha)*i - sin(alpha)*q
   q  = sin(alpha)*i + cos(alpha)*q
   i  = i1

   if defined(ires) eq 0 then begin
      n = n_elements(i)
      ires = dblarr( nflux, n)
      qres = dblarr( nflux, n)
      xcres = dblarr( nflux)
      ycres = dblarr( nflux)
      rres  = dblarr( nflux)
   endif
   ires[ iflux,*] = i
   qres[ iflux,*] = q
   xcres[iflux,*] = xc
   ycres[iflux,*] = yc
   rres[ iflux,*] = r
endfor

;; Define the ouput LF time vector:
;; need to add 0.5/f_sampling to compute the time corresponding to the average
;; of the 20-40 points interval taken in iqdidq to produce a final sample from
;; the "oversampled" i,q timelines.
ni = n_elements(toi_rf_planet_hz)
time_lf = dindgen(ni)/f_sampling + 0.5/f_sampling

;; Circle plot
nphi = 10000
phi = dindgen(nphi)/nphi*2*!dpi
cosphi = cos(phi)
sinphi = sin(phi)
j = dcomplex(0,1)

xra = [-0.5,10]
yra = [-0.5,5]

position = [0.1, 0.1, 0.95, 0.95]
position1 = [position[0], 0.83, 0.95, 0.95]
wind, 1, 1, /free, /large
outplot, file='circle_zres', png=png, ps=ps
plot, ires[1,*], qres[1,*], /xs, /ys, $
      xtitle='In phase', ytitle='In Quadrature', psym=1, /iso, $
      xra=xra, yra=yra, /nodata, position=position
oplot, [0,0], [-1,1]*1d10, line=2
oplot, [-1,1]*1d10, [0,0], line=2
oplot, [1,1], [-1,1]*1d10, line=2
for iflux=nflux-1, 0, -1 do begin
   oplot, ires[iflux,*], qres[iflux,*], psym=1, col=flux_col[iflux]
   oplot, xcres[iflux] + rres[iflux]*cosphi, ycres[iflux]+rres[iflux]*sinphi, col=flux_col[iflux]

   alpha = atan(ycres[iflux],xcres[iflux])
   Inorm = -1/(2.*rres[iflux]) * ((ires[iflux,*]-xcres[iflux])*cos(alpha) + $
                                  (Qres[iflux,*]-ycres[iflux])*sin(alpha)) + 1/2.
   Qnorm = -1/(2.*rres[iflux]) * (-(ires[iflux,*]-xcres[iflux])*sin(alpha) + $
                                  (Qres[iflux,*]-ycres[iflux])*cos(alpha))
   Znorm = Inorm + j*Qnorm
   Zres  = 1d0/Znorm

   oplot, float(Zres), imaginary(Zres), psym=8, syms=0.5, col=flux_col[iflux]
endfor
legendastro, strtrim( in_flux_list,2)+" Jy", col=flux_col, line=0

iflux=nflux-1
xra = [14,16]
dd = sqrt( (ires[iflux,*]-xcres[iflux])^2 + $
           (qres[iflux,*]-ycres[iflux])^2) - rres[iflux]
plot, time_lf, dd/rres[iflux], /xs, xra=xra, $
      position=position1, /noerase, /ys
for iflux=nflux-1, 0, -1 do begin
   dd = sqrt( (ires[iflux,*]-xcres[iflux])^2 + $
                (qres[iflux,*]-ycres[iflux])^2) - rres[iflux]
   oplot, time_lf, dd/rres[iflux], col=flux_col[iflux]
endfor
outplot, /close, /verb

end
