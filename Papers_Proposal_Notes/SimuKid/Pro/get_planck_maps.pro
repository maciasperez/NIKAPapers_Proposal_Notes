
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

;; Dust at 353GHz
planck_dust = mrdfits("$SK_DIR/Maps/COM_CompMap_DustPol-commander_1024_R2.00.fits", 1, hdust)
print, sxpar(hdust,"NU_REF")
lambda = !const.c/353.d9 * 1d6 ; microns
planck_dust.q_ml_full /= 1000.d0 ; microK_RJ to mK_RJ
planck_dust.u_ml_full /= 1000.d0
convert_millik_megajy, lambda, planck_dust.q_ml_full, q_megajy, /rj
convert_millik_megajy, lambda, planck_dust.u_ml_full, u_megajy, /rj
convert_megajy_millik, lambda, q_megajy, q_dust_ref, /cmb_th
convert_megajy_millik, lambda, u_megajy, u_dust_ref, /cmb_th
q_dust_ref *= 1000 ; mK to uK
u_dust_ref *= 1000 ; mK to uK
ud_grade, q_dust_ref, q_dust_ref_out, nside_out=nside, order_in='nested', order_out='ring'
ud_grade, u_dust_ref, u_dust_ref_out, nside_out=nside, order_in='nested', order_out='ring'
q_dust_ref = q_dust_ref_out
u_dust_ref = u_dust_ref_out

planck_t_dust = mrdfits("$SK_DIR/Maps/COM_CompMap_ThermalDust-commander_2048_R2.00.fits", 1, htdust)
print, sxpar( htdust, "NU_REF")
lambda = !const.c/545.d9 * 1d6 ; microns
mollview, planck_t_dust.i_ml_full, /nest
ud_grade, planck_t_dust.i_ml_full, i_dust_ref, nside_out=nside, order_in='nested', order_out='ring'
convert_millik_megajy, lambda, i_dust_ref, i_megajy, /rj
convert_megajy_millik, lambda, i_megajy, i_dust_ref, /cmb_th
i_dust_ref *= 1000 ; mK to uK

;; Synchrotron
planck_sync = mrdfits("$SK_DIR/Maps/COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits", 1, hsync)
print, sxpar( hsync, "NU_REF")
lambda = !const.c/30.d9 * 1d6 ; microns

;; uK_RJ to mK_RJ
planck_sync.q_ml_full /= 1000.d0 ; uK to mK
planck_sync.u_ml_full /= 1000.d0
convert_millik_megajy, lambda, planck_sync.q_ml_full, q_megajy, /rj
convert_megajy_millik, lambda, q_megajy, q_sync_ref, /cmb_th
convert_millik_megajy, lambda, planck_sync.u_ml_full, u_megajy, /rj
convert_megajy_millik, lambda, u_megajy, u_sync_ref, /cmb_th
q_sync_ref *= 1000 ; mK to uK
u_sync_ref *= 1000 ; mK to uK
ud_grade, q_sync_ref, q_sync_ref_out, nside_out=nside, order_in='nested', order_out='ring'
ud_grade, u_sync_ref, u_sync_ref_out, nside_out=nside, order_in='nested', order_out='ring'
q_sync_ref = q_sync_ref_out
u_sync_ref = u_sync_ref_out

planck_t_sync = mrdfits("$SK_DIR/Maps/COM_CompMap_Synchrotron-commander_0256_R2.00.fits", 1, htsync)
print, sxpar( htsync, "NU_REF")
lambda = !const.c/408.d6 * 1d6 ; microns
ud_grade, planck_t_sync.i_ml, i_sync_ref, nside_out=nside, order_in='nested', order_out='ring'
i_sync_ref /= 1000.d0 ; uK to mK
convert_millik_megajy, lambda, i_sync_ref, i_megajy, /rj
convert_megajy_millik, lambda, i_megajy, i_sync_ref, /cmb_th
i_sync_ref *= 1000 ; mK to uK

if not keyword_set(noplot) then begin
   mollview, i_dust_ref, title='I dust (from 545 GHz)'
   mollview, q_dust_ref, title='Q dust (from 353 GHz)'
   mollview, u_dust_ref, title='U dust (from 353 GHz)'
   mollview, i_sync_ref, title='I sync_ref (from 408 MHz)'
   mollview, q_sync_ref, title='Q sync_ref (from 30 GHz)'
   mollview, u_sync_ref, title='U sync_ref (from 30 GHz)'
endif

end
