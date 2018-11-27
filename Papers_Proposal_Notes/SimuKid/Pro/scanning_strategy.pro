;; from plot_galaxy.pro
      ;; from plot_planet_2_test.pro
          ;; Scan strat (from core_scan.F90)
;;----------------------------------------------

col_rf = 70
col_cf = 250

png = 0 ;1
ps  = 0 ;1

;; Scanning strategy parameters
;;Planck
manip = 'Planck'
alpha_deg = 7.d0
beta_deg  = 85.d0
gamma_deg = 0.d0
t_pre_sec = 365./2*86400.d0
;; ;; t_spin_sec = 60.

scan_strat, alpha_deg, beta_deg, gamma_deg, t_pre_sec, $
            toi_dust_jy=toi_dust_planck_jy, toi_dipole_jy=toi_dipole_planck_jy, $
            toi_dust_dipole_jy=toi_dust_dipole_planck_jy

;; EPIC
manip = 'EPIC'                  ;'Pol. Sat.'
alpha_deg   = 55.d0
beta_deg    = 45.d0
gamma_deg   = 0.d0
t_pre_sec = 2*3600.d0           ; sec ; precession around Sun-Earth
;; ;; t_spin_sec  = 301.d0 ; 61.d0 ; sec ; rotation about spin axis

scan_strat, alpha_deg, beta_deg, gamma_deg, t_pre_sec, $
            toi_dust_jy=toi_dust_epic_jy, toi_dipole_jy=toi_dipole_epic_jy, $
            toi_dust_dipole_jy=toi_dust_dipole_epic_jy


my_multiplot, 1, 2, pp, pp1, /rev

if ps eq 0 then wind, 1, 1, /free, /large
outplot, file='histo_galaxy_dipole', png=png, ps=ps
np_histo, toi_dust_dipole_planck_jy, xtitle='Flux (Jy)', position=pp1[0,*]
np_histo, toi_dust_dipole_epic_jy, xtitle='Flux (Jy)', position=pp1[1,*], /noerase
outplot, /close, /verb




end
