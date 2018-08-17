
pro get_planck_maps, nside, $
                     i_cmb, q_cmb, u_cmb, $
                     i_dust_ref, q_dust_ref, u_dust_ref, $
                     i_sync_ref, q_sync_ref, u_sync_ref, noplot=noplot

map_dir = '$SK_DIR/Maps'

planck_cmb = mrdfits("$SK_DIR/Maps/COM_CMB_IQU-commander_1024_R2.02_full.fits", 1, hcmb)
;; K to microK
i_cmb = planck_cmb.i_stokes * 1d6
q_cmb = planck_cmb.q_stokes * 1d6
u_cmb = planck_cmb.u_stokes * 1d6
ud_grade, i_cmb, i_cmb, nside_out=nside, order_in='nested', order_out='ring'
ud_grade, q_cmb, q_cmb, nside_out=nside, order_in='nested', order_out='ring'
ud_grade, u_cmb, u_cmb, nside_out=nside, order_in='nested', order_out='ring'

;; Dust Q and U at 353GHz
planck_dust = mrdfits("$SK_DIR/Maps/COM_CompMap_DustPol-commander_1024_R2.00.fits", 1, hdust)
print, sxpar(hdust,"NU_REF")
print, sxpar(hdust, "TUNIT1")
ud_grade, planck_dust.q_ml_full, q_dust_ref, nside_out=nside, order_in='nested', order_out='ring'
ud_grade, planck_dust.u_ml_full, u_dust_ref, nside_out=nside, order_in='nested', order_out='ring'
;;lambda = !const.c/353.d9 * 1d6 ; microns
;;q_dust_ref /= 1000.d0 ; microK_RJ to mK_RJ
;;u_dust_ref /= 1000.d0
;;;; Convert from RJ to thermo
;;convert_millik_megajy, lambda, q_dust_ref, q_megajy, /rj
;;convert_millik_megajy, lambda, u_dust_ref, u_megajy, /rj
;;convert_megajy_millik, lambda, q_megajy, q_dust_ref, /cmb_th
;;convert_megajy_millik, lambda, u_megajy, u_dust_ref, /cmb_th
;;q_dust_ref *= 1000              ; mK to uK
;;u_dust_ref *= 1000              ; mK to uK

;; Dust T at 545GHz
planck_t_dust = mrdfits("$SK_DIR/Maps/COM_CompMap_ThermalDust-commander_2048_R2.00.fits", 1, htdust)
print, sxpar( htdust, "NU_REF")
print, sxpar( htdust, "TUNIT1")
ud_grade, planck_t_dust.i_ml_full, i_dust_ref, nside_out=nside, order_in='nested', order_out='ring'
;;lambda = !const.c/545.d9 * 1d6 ; microns
;;i_dust_ref /= 1000 ; uK to mK
;;convert_millik_megajy, lambda, i_dust_ref, i_megajy, /rj
;;convert_megajy_millik, lambda, i_megajy, i_dust_ref, /cmb_th
;;i_dust_ref *= 1000              ; mK to uK

;; Synchrotron
planck_sync = mrdfits("$SK_DIR/Maps/COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits", 1, hsync)
print, sxpar( hsync, "NU_REF")
print, sxpar( hsync, "TUNIT1")
ud_grade, planck_sync.q_ml_full, q_sync_ref, nside_out=nside, order_in='nested', order_out='ring'
ud_grade, planck_sync.u_ml_full, u_sync_ref, nside_out=nside, order_in='nested', order_out='ring'
;;lambda = !const.c/30.d9 * 1d6 ; microns
;;q_sync_ref /= 1000.d0 ; uK to mK
;;u_sync_ref /= 1000.d0
;;convert_millik_megajy, lambda, q_sync_ref, q_megajy, /rj
;;convert_millik_megajy, lambda, u_sync_ref, u_megajy, /rj
;;convert_megajy_millik, lambda, q_megajy, q_sync_ref, /cmb_th
;;convert_megajy_millik, lambda, u_megajy, u_sync_ref, /cmb_th
;;q_sync_ref *= 1000 ; mK to uK
;;u_sync_ref *= 1000 ; mK to uK

planck_t_sync = mrdfits("$SK_DIR/Maps/COM_CompMap_Synchrotron-commander_0256_R2.00.fits", 1, htsync)
print, sxpar( htsync, "NU_REF")
print, sxpar( htsync, "TUNIT1")
ud_grade, planck_t_sync.i_ml, i_sync_ref, nside_out=nside, order_in='nested', order_out='ring'
;;lambda = !const.c/408.d6 * 1d6 ; microns
;;i_sync_ref /= 1000.d0 ; uK to mK
;;convert_millik_megajy, lambda, i_sync_ref, i_megajy, /rj
;;convert_megajy_millik, lambda, i_megajy, i_sync_ref, /cmb_th
;;i_sync_ref *= 1000 ; mK to uK

if not keyword_set(noplot) then begin
   mollview, i_dust_ref, title='I dust (from 545 GHz)'
   mollview, q_dust_ref, title='Q dust (from 353 GHz)'
   mollview, u_dust_ref, title='U dust (from 353 GHz)'
   mollview, i_sync_ref, title='I sync_ref (from 408 MHz)'
   mollview, q_sync_ref, title='Q sync_ref (from 30 GHz)'
   mollview, u_sync_ref, title='U sync_ref (from 30 GHz)'
endif

end
