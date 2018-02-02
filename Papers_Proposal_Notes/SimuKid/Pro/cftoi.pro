
function cftoi, i, q, di, dq, delta_f, xc, yc, r, polydeg=polydeg, stop=stop


if not keyword_set(polydeg) then polydeg = 2

j = dcomplex(0,1)

alpha = atan(yc,xc)

;; ***** Changing the sign, NP, Jan. 18th ****
Inorm = -1/(2.*r) * ( (I-xc)*cos(alpha) + (Q-yc)*sin(alpha)) + 1/2.
Qnorm = -1/(2.*r) * (-(I-xc)*sin(alpha) + (Q-yc)*cos(alpha))
Znorm = Inorm + j*Qnorm
;circlefit, inorm, qnorm, xc2, yc2, r2
;print, "xc2, yc2, r2: ", xc2, yc2, r2
;stop

;; ***** Changing the sign, NP, Jan. 18th ****
dInorm = -1/(2.*r)*( dI*cos(alpha) + dQ*sin(alpha))
dQnorm = -1/(2.*r)*(-dI*sin(alpha) + dQ*cos(alpha))
dZnorm = dInorm +j*dQnorm

Zres  = 1d0/Znorm
dZres = -dZnorm/(Znorm^2)

dy3 = imaginary(dZres)
y3  = imaginary(Zres)
y3r = double(Zres)

deriv = delta_f/dy3
;; stop
;; noise and fit of the transformed circle 
m = yc/xc
xres = ((xc+m*yc)-r*sqrt(1+m^2))/(1+m^2)
yres = m*xres
xop  = ((xc+m*yc)+r*sqrt(1+m^2))/(1+m^2)
yop  = m*xop

dnoise = 1/sqrt((i-xop)^2 + (q-yop)^2)
noisefit = dnoise/dy3^2

Rn  =  poly_fit( y3, deriv, polydeg, /double, yfit=yfitRn, measure_errors=noisefit, status=status)

if keyword_set(stop) then begin
   wind, 1, 1, /free, /large
   !p.multi=[0,2,2]
   plot, y3, deriv,psym=4, /xs, /ys
   oplot, y3, deriv, psym=4, col=250
   oplot, y3, yfitRn, col=150
   legendastro, 'Where is yfitRn ????!!!'
   plot, y3, yfitrn, /xs, /ys
   plot, y3, deriv, /xs, /ys
   !p.multi=0
   stop
endif

;; Integrate the fitted polynomial
pn = i*0.d0
for ii=0, polydeg do pn += Rn[ii]/(ii+1.d0)*y3^(ii+1)

return, reform(pn)
end
