
;; Xavier's estimation of Planet fluxes: take this as an order of
;; magnitude cross checks of our unit conversions
;;
;; see Labtools/FXD/Run9/Misc.scr
;;-----------------------------------------------

;; Typical Uranus
T_Uranus = 93  ; K
lambda = [1.2, 2]*1d-3 ; NIKA
size_Uranus = sqrt(3.68* 3.57) ; Oct 2014
nubnu_planck, lambda, T_Uranus, nubnu
conv_nubnu_bnu, lambda, nubnu, bnu
flux_source = bnu * 1E6 * !pi/4.* (size_Uranus/3600./!radeg)^2 ;
print, ';; Uranus flux (NIKA bands) (Jy): ', flux_source
;; Uranus flux (NIKA bands) (Jy):        40.520387       14.974997

;; CMB Dipole in a typical 7arcmin beam
print, ""
T_dipole = 2.726 + [-1,1]*3d-3/2.
lambda = 3d-3 ; <-- 100GHz
fwhm = 7. ; arcmin
beam_size = sqrt(8.*(fwhm*!fwhm2sigma*!arcmin2rad)^2) ; effective diameter
nubnu_planck, lambda, T_dipole, nubnu
conv_nubnu_bnu, lambda, nubnu, bnu
flux_dipole = bnu * 1E6 * !pi/4.* beam_size^2 ;
print, ';; CMB Dipole flux variation at 100GHz in a 7 arcmin beam (Jy): ', max(flux_dipole)-min(flux_dipole)
;; CMB Dipole flux variation at 100GHz in a 7 arcmin beam (Jy):        3.3621599

;; Galactic Dust at 353GHz
planck_dust = mrdfits("$SK_DIR/Maps/COM_CompMap_DustPol-commander_1024_R2.00.fits", 1, hdust)
print, sxpar(hdust,"NU_REF")
print, sxpar(hdust,"TUNIT1")
lambda_353 = !const.c/353.d9 * 1d6 ; 353
;; convert from uK_RJ to mK_RJ and from mK_RJ to mK_CMB
convert_millik_megajy, lambda_353, planck_dust.q_ml_full/1000., q_megajy, /rj
convert_millik_megajy, lambda_353, planck_dust.u_ml_full/1000., u_megajy, /rj
convert_megajy_millik, lambda_353, q_megajy, q_dust_ref, /cmb_th
convert_megajy_millik, lambda_353, u_megajy, u_dust_ref, /cmb_th
;; Nest to ring
ud_grade, q_dust_ref, q_dust_ref, nside_out=nside, order_in='nested', order_out='ring'
ud_grade, u_dust_ref, u_dust_ref, nside_out=nside, order_in='nested', order_out='ring'

planck_t_dust = mrdfits("$SK_DIR/Maps/COM_CompMap_ThermalDust-commander_2048_R2.00.fits", 1, htdust)
print, sxpar( htdust, "NU_REF")
lambda_545 = !const.c/545.d9 * 1d6 ; microns
ud_grade, planck_t_dust.i_ml_full, i_dust_ref, nside_out=nside, order_in='nested', order_out='ring'
;; convert from uK_RJ to mK_RJ and from mK_RJ to mK_CMB
convert_millik_megajy, lambda_545, i_dust_ref/1000., i_megajy, /rj
convert_megajy_millik, lambda_545, i_megajy, i_dust_ref, /cmb_th

;; Extrapolate Idust(545GHz) to 353GHz
beta_dust = 1.5d0
i_dust_ref *= (353./545)^beta_dust

i_dust_ref /= 1000.d0 ; mK to K
q_dust_ref /= 1000.d0 ; mK to K
u_dust_ref /= 1000.d0 ; mK to K

print, ";; minmax(i_dust_ref) (K): ", minmax(i_dust_ref)
print, ";; minmax(q_dust_ref) (K): ", minmax(q_dust_ref)
print, ";; minmax(u_dust_ref) (K): ", minmax(u_dust_ref)
;; minmax(i_dust_ref) (K):    3.7330228e-05       23.881363
;; minmax(q_dust_ref) (K):    -0.0053505400    0.0062592595
;; minmax(u_dust_ref) (K):    -0.0044205532     0.016788402
;stop

;; Apply Uranus type flux derivation to our diffuse emissions
temp_name = ['Qdust', 'Udust', 'Q(10%) pol']
t_dust = 10. ; typical...?
dT = [ max(abs(q_dust_ref)), max(abs(u_dust_ref)), 0.1*t_dust]/2.
lambda = 3d-3 ; <-- 100GHz
fwhm = 7. ; arcmin
;; take effective disk diameter d such that pi*d^2/4 = 2*!dpi*sigma^2
beam_size = sqrt( 8.*(fwhm*!fwhm2sigma*!arcmin2rad)^2)
for i=0, n_elements(dt)-1 do begin
   T = t_dust + dt[i]
   nubnu_planck, lambda, T, nubnu
   conv_nubnu_bnu, lambda, nubnu, bnu
   flux1 = bnu * 1E6 * !pi/4.* beam_size^2 ;

   T = t_dust - dt[i]
   nubnu_planck, lambda, T, nubnu
   conv_nubnu_bnu, lambda, nubnu, bnu
   flux2 = bnu * 1E6 * !pi/4.* beam_size^2

   print, ';; '+temp_name[i]+" flux variation in a 7 arcmin beam (Jy): "+strtrim(abs(flux2-flux1),2)
endfor
;; Qdust flux variation in a 7 arcmin beam (Jy): 8.8507145
;; Udust flux variation in a 7 arcmin beam (Jy): 23.739127
;; Q(10%) pol flux variation in a 7 arcmin beam (Jy): 1413.9528


;; Check on the converted maps
print, ";; minmax(q_megajy)*2.d0*!dpi*(fwhm*!fwhm2sigma*!arcmin2rad)^2: "
print, ";; ", minmax(q_megajy)*2.d0*!dpi*(fwhm*!fwhm2sigma*!arcmin2rad)^2*1d6
print, ";; minmax(u_megajy)*2.d0*!dpi*(fwhm*!fwhm2sigma*!arcmin2rad)^2: "
print, ";; ", minmax(u_megajy)*2.d0*!dpi*(fwhm*!fwhm2sigma*!arcmin2rad)^2*1d6
;; minmax(q_megajy)*2.d0*!dpi*(fwhm*!fwhm2sigma*!arcmin2rad)^2: 
;;       -7.4628997       8.7303758
;; minmax(u_megajy)*2.d0*!dpi*(fwhm*!fwhm2sigma*!arcmin2rad)^2: 
;;       -6.1657598       23.416358

;; Dipole flux:
lambda_microns = lambda*1.d6
convert_millik_megajy, lambda_microns, 1000*(2.726 + 3d-3/2.), i_dipole_megajy1, /cmb
convert_millik_megajy, lambda_microns, 1000*(2.726 - 3d-3/2.), i_dipole_megajy2, /cmb
print, ";; Dipole flux in a 7 arcmin beam (Jy): ", 1.d6*(i_dipole_megajy1-i_dipole_megajy2)*2.d0*!dpi*(fwhm*!fwhm2sigma*!arcmin2rad)^2
;; Dipole flux in a 7 arcmin beam (Jy):    3.3623570

;; Analytical conversion
t_cmb = 2.726
freq = !const.c/lambda
print, "freq: ", freq
;; W.m^-2.sr^-1.Hz^-1 /K
dBBdT = 2*!const.h^2*freq^4/!const.k/!const.c^2/t_cmb^2 * $
        EXP(!const.h/!const.k*freq/t_cmb) / (EXP(!const.h/!const.k*freq/t_cmb) - 1)^2

;; (W.m^-2.Hz^-1/sr)/K to (Jy/sr)/K (cross-checked at 100GHz with planck_g() in convert_millik_megajy)
dBBdT = dBBdT/1d-26

flux1 = dbbdt * (3+3.d-3/2.)
flux2 = dbbdt * (3-3.d-3/2.)
print, ";; Dipole amplitude flux in 7 arcmin beam (analytical) (Jy): ", $
       (flux1-flux2)*2*!dpi*(fwhm*!fwhm2sigma*!arcmin2rad)^2
;; Dipole amplitude flux in 7 arcmin beam (analytical) (Jy):        3.3623260


end
