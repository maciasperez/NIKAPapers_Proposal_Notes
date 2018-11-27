pro scan_strat, alpha_deg, beta_deg, gamma_deg, t_pre_sec, t_spin_sec, $
                toi_dust_jy=toi_dust_jy, toi_dipole_jy=toi_dipole_jy, toi_dust_dipole_jy=toi_dust_dipole_jy



nside = 2048
;nside_out = 2048
fwhm_arcmin = 5. ; arcmin
n_harmonics = 8 ;6  ; hwp template harmonics

;; Number of harmonics of the HWP that we want to determine (implies f_sampling
;; >= 2*n_harmonics*hwp_rot_freq)
;; We first fix the HWP rotation frequency
hwp_rot_freq = 3.d0 ; 3.d0 like NIKA
f_sampling = long( (12.*hwp_rot_freq) > (2*n_harmonics*hwp_rot_freq))
npts_per_fwhm_min = 3
vmax = fwhm_arcmin/npts_per_fwhm_min * 4 * hwp_rot_freq ; arcmin/s

;; To produce signal timelines from pixelized maps
f_nyquist = 1./(fwhm_arcmin/3./vmax)



t_spin_sec = sin(beta_deg*!dtor)/(vmax/(360.*60.)-sin(alpha_deg*!dtor)/t_pre_sec)
print, "t_spin_sec: ", t_spin_sec

;; Scan parameters and their units
alpha = alpha_deg * !dtor
beta  = beta_deg  * !dtor
gamma = gamma_deg * !dtor
;; print, "T_alpha (s): ", t_pre_sec
;; print, "T_beta (s): ", t_spin_sec
;; print, "npts_per_fwhm_min: ", npts_per_fwhm_min
;; print, "vmax (arcmin/s): ", vmax
;; print, "f_sampling (for polar): ", f_sampling
;; print, "f_nyquist (to scan the input maps): ", f_nyquist

;; Oversample to prepare I, Q, dI, dQ
fsampl_factor = 20L
f_modulation = fsampl_factor*long(f_sampling)

obs_time_hour = 10.d0/60.d0 ; 1./60d0
nsn = round(obs_time_hour*3600.d0*f_modulation)
time_modulation = dindgen(nsn)/f_modulation
nsn_nyquist = round( obs_time_hour*3600.d0*f_nyquist)

;; Retrieve sky maps
if defined(dust_mjysr) eq 0 then init_maps, nside, dust_mjysr, dipole_mjysr

;; mollview, dipole_mjysr, title='Dipole (Mjy/sr)'
;; mollview, dust_mjysr, title='map (Mjy/sr)'
;; mollview, dipole_mjysr+dust_mjysr, title='Raw + dipole (Mjy/sr)'

;; Pointing
get_pointing, nside, t_pre_sec, t_spin_sec, beta, nsn_nyquist, f_nyquist, alpha, ipix_vec

;; Timelines of Galaxy (dust) and dipole in Jy at f_nyquist
time_nyquist = dindgen(nsn_nyquist)/f_nyquist
omega_beam_sr    = 2.d0*!dpi*(fwhm_arcmin*!fwhm2sigma*!arcmin2rad)^2
toi_dust_jy      = omega_beam_sr * dust_mjysr[  ipix_vec] * 1.d6
toi_dipole_jy    = omega_beam_sr * dipole_mjysr[ipix_vec] * 1.d6
toi_dust_dipole_jy       = omega_beam_sr * (dust_mjysr + dipole_mjysr)[ipix_vec] * 1.d6

;; Perform interpolation at f_modulation
toi_dust_jy   = interpol( toi_dust_jy,   time_nyquist, time_modulation, /spline)
toi_dipole_jy = interpol( toi_dipole_jy, time_nyquist, time_modulation, /spline)

;; mask[ipix_mask] = 1
;; mollview, mask, max=2, title=manip


toi_dust_dipole_jy = toi_dust_jy + toi_dipole_jy


end
