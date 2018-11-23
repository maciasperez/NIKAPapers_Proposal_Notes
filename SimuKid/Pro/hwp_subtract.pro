; from nk_hwp_rm4.pro
;----------------------

pro hwp_subtract, toi, f_sampling, omega_rad, n_harmonics, toi_out, hwp_beta_out, toi_mask=toi_mask, $
                  drift=drift

nsn = n_elements(toi)
if not keyword_set(toi_mask) then begin
   wkeep = lindgen(nsn)
   nwmask = 0
endif else begin
   wkeep = where( toi_mask eq 1, nwkeep, compl=wmask, ncompl=nwmask)
   if nwkeep eq 0 then begin
      message, /info, "No pixel far from the source"
      stop
   endif
endelse

;; Mimic the interpolation on NIKA points but here the acquisition is not
;; synchronous with the HWP rotation
womega = where( abs(omega_rad) le abs(omega_rad[1]-omega_rad[0]), nwomega)
;temp_toi = toi[0:womega[nwomega-1]]

index = lindgen(nsn)

;; Subtract HWPSS and decorrelate atmosphere
;; my_interpol = dblarr(nsn)
;; for i=0, nwomega-2 do begin
;;    dx    = double(womega[i+1] - womega[i])
;;    dy    = toi[womega[i+1]] - toi[womega[i]]
;;    n     = womega[i+1]-womega[i]+1
;;    xx    = dindgen(n)
;;    my_interpol[womega[i]:womega[i+1]] = toi[womega[i]] + (dy/dx)*xx
;; endfor

my_interpol = interpol( toi[womega], womega, index)
;; wind, 1, 1, /free, /xlarge
;; plot, toi, xra=[0,1000],/xs, /ys
;; oplot, womega, toi[womega], psym=4, col=70
;; oplot, my_interpol, col=250
;; stop

;; nk_hwp_rm_4_sub:
ncoeff = 2 + 4*n_harmonics
time = dindgen(nsn)/f_sampling

;; Build model harmonics

amplitudes = dblarr( ncoeff-2)
outtemp = dblarr( ncoeff-2, nsn)
for i=0, n_harmonics-1 do begin
   outtemp[ i*4,     *] =      cos( (i+1)*omega_rad)
   outtemp[ i*4 + 1, *] = time*cos( (i+1)*omega_rad)
   outtemp[ i*4 + 2, *] =      sin( (i+1)*omega_rad)
   outtemp[ i*4 + 3, *] = time*sin( (i+1)*omega_rad)
endfor


if drift eq 1 then begin 
   amplitudes = dblarr( ncoeff)
   outtemp = dblarr( ncoeff, nsn)
   outtemp[0,*] = 1.0d0
   outtemp[1,*] = time
   for i=0, n_harmonics-1 do begin
      outtemp[ 2 + i*4,     *] =   cos( (i+1)*omega_rad)
      outtemp[ 2 + i*4 + 1, *] = time*cos( (i+1)*omega_rad)
      outtemp[ 2 + i*4 + 2, *] =   sin( (i+1)*omega_rad)
      outtemp[ 2 + i*4 + 3, *] = time*sin( (i+1)*omega_rad)
   endfor
endif




;; Subtract my_interpol and the constant to prepare HWP
;; fit outside the source
;; The constant is computed on the correct range but
;; can be subtracted everywhere
toi1 = toi - my_interpol
if nwmask ne 0 then toi1[wmask] = !values.d_nan
a = avg( toi1, /nan)
toi1 -= a
if nwmask ne 0 then toi1[wmask] = 0.d0 ; remove the NaN for the next loop

;; Fit the HWPSS outside the source and on valid samples
temp = outtemp                  ; init
if nwmask ne 0 then temp[*,wmask] = 0.d0
ata   = matrix_multiply( temp, temp, /btranspose)
atam1 = invert(ata)
atd        = transpose(temp) ## toi1
amplitudes = atam1##atd
hwp_beta_out    = outtemp##amplitudes

;stop

toi_out = toi - hwp_beta_out - a


end
