
pro freso2toi, freso_in, kid_model, delta_f, toi_rf, toi_cf, epsilon_qi=epsilon_qi, $
               i=i, q=q, di=di, dq=dq, npts_avg=npts_avg, test=test, apply_shift=apply_shift, $
               stop=stop, xc=xc, yc=yc, r=r
  
f_m = kid_model.f0 - delta_f/2.
f_p = kid_model.f0 + delta_f/2.

freso = freso_in + kid_model.f0

s21_m = resonance( kid_model, freso, f_m, epsilon_qi=epsilon_qi)
s21_p = resonance( kid_model, freso, f_p, epsilon_qi=epsilon_qi)

nsn = n_elements(freso)
pt = lindgen(nsn)
s21 = dcomplexarr(2*nsn)
s21[2*pt]   = s21_p
s21[2*pt+1] = s21_m

iqdidq, kid_model, s21, I, Q, dI, dQ, npts_avg=npts_avg

;; method rf
toi_rf = -rf(I, Q, dI, dQ, delta_f, apply_shift=apply_shift)
npn = n_elements(toi_rf)
if not keyword_set(test) then toi_rf = toi_rf - toi_rf[npn-1]

;; method cf
; iqfit = linfit(i,q, chisqr=khi2)
circlefit, I, Q, xc, yc, r, avgdist, zzmin=1d-6
toi_cf = -cftoi(i, q, di, dq, delta_f, xc, yc, r, polydeg = 2, stop=stop)
if not keyword_set(test) then toi_cf = toi_cf - toi_cf[npn-1]

;; ;;===============
;; nphi = 100000L
;; phi = dindgen(nphi)/(nphi-1)*2*!dpi
;; wind, 1, 1, /free
;; plot, i, q, /xs, /ys, psym=1, xtitle='I', ytitle='Q'
;; oplot, xc + r*cos(phi), yc + r*sin(phi), col=250;, psym=8
;; legendastro, ["xc: "+strtrim(xc,2), $
;;               "yc: "+strtrim(yc,2), $
;;               "r: "+strtrim(r,2)]
;; stop


end

