
pro get_pointing, nside, t_pre_sec, t_spin_sec, beta, nsn_nyquist, f_nyquist, alpha, ipix_vec

;; Precession about the sun-spacecraft axis
omega_alpha_pre = 2.d0*!dpi/t_pre_sec

; Spin about the symetry axis of the payload
omega_beta_spin = 2.d0*!dpi/t_spin_sec

; earth around the sun
omega_earth = 2.d0*!dpi/ 60.d0 / 60.d0 / 24.d0 / 365.d0

v0 = dblarr(3)
;;direction de visee du detecteur a t=0
V0[0] = cos(beta)
V0[1] = 0.d0
V0[2] = sin(beta)

;; ipix_vec = lonarr(nsn)
;; ipix_mask = lonarr(nsn)

ipix_vec  = lonarr(nsn_nyquist)
ipix_mask = lonarr(nsn_nyquist)

nside_mask = 128
mask = fltarr( nside2npix(nside_mask))

;; for isn=0L, nsn-1 do begin
;;    t   = float(isn)/f_modulation

for isn=0L, nsn_nyquist-1 do begin
   t = float(isn)/f_nyquist

   phi_earth = (omega_earth     * t) mod (2.d0*!dpi) ; earth to sun angle
   beta_spin = (omega_beta_spin * t) mod (2.d0*!dpi) ; spin about the symetry axis of the sat.
   alpha_pre = (omega_alpha_pre * t) mod (2.d0*!dpi) ; precession about the sun-spacecraft axis

   ;; On place l'axe de spin dans le repere heliocentrique
   spin_axis = [cos(alpha), 0.d0, sin(alpha)]

   ;; le detecteur fait l'angle beta par rapport a l'axe de spin
   ;; => dans la base (spin_axis, v, w):
   v = [-spin_axis[2], 0., spin_axis[0]]
   w = crossp( spin_axis, v)
   
   ;; Matrice de Passage P
   P = [[spin_axis[0], v[0], w[0]], $
        [spin_axis[1], v[1], w[1]], $
        [spin_axis[2], v[2], w[2]]]

   ;; On fait tourner le detecteur autour de l'axe de spin du satellite
   det = [cos(beta), 0.d0, sin(beta)]
   r = [[1.d0, 0.d0, 0.d0], $
        [0.d0, cos(beta_spin), -sin(beta_spin)], $
        [0.d0, sin(beta_spin),  cos(beta_spin)]]
   d = r##det

   ;; On exprime d dans la base heliocentrique
   d = P##d

   ;; On fait precesser le sat dans cette base autour de x = Sun-Earth
   R = [[1.d0, 0.d0, 0.d0], $
        [0.d0, cos(alpha_pre), -sin(alpha_pre)], $
        [0.d0, sin(alpha_pre),  cos(alpha_pre)]]
   d = R##d

   ;; Rotation d'axe Z pour decrire la rotation de la terre autour du soleil
   R = [[cos(phi_earth), -sin(phi_earth), 0.d0], $
        [sin(phi_earth),  cos(phi_earth), 0.d0], $
        [0.d0, 0.d0, 1.d0]]
   d = R##d
   
   vec2pix_ring, nside, d, ipix
   ipix_vec[isn] = ipix

   vec2pix_ring, nside_mask, d, ipix
   ipix_mask[isn] = ipix
endfor

;mask[ipix_mask] = 1
;mollview, mask, max=2, title=manip

end
