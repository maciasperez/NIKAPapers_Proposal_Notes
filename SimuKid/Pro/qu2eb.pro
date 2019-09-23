
pro qu2eb, map_q, map_u, res_arcmin, map_e, map_b

s = size(map_q)
n = s[1]
lmap_rad = float(n*res_arcmin*!arcmin2rad)

amn_q = fft( map_q, /double)
amn_u = fft( map_u, /double)
amn_e = amn_q*0.d0
amn_b = amn_q*0.d0

ic = complex( 0.00d0, 1.00d0)

for im=0L, n-1 do begin
   if im le n/2 then m1 = float(im) else m1 = float(n-im)
   for in=0L, n-1 do begin
      if in le n/2 then n1 = float(in) else n1 = float(n-in)
      um = m1/float(Lmap_rad)
      un = n1/float(Lmap_rad)
      u2 = um^2 + un^2
      if u2 eq 0 then begin
         amn_e[im,in] = 0.d0
         amn_b[im,in] = 0.d0
      endif else begin
         ;amn_q[im,in] = amn_e[im,in]*(um^2 - un^2)/u2 - amn_b[im,in]*2.0d0*um*un/u2
         ;amn_u[im,in] = amn_e[im,in]*2.0d0*um*un/u2   + amn_b[im,in]*(um^2 - un^2)/u2
         amn_e[im,in] =  amn_q[im,in]*(um^2 - un^2)/u2 + amn_u[im,in]*2.0d0*um*un/u2
         amn_b[im,in] = -amn_q[im,in]*2.0d0*um*un/u2   + amn_u[im,in]*(um^2 - un^2)/u2
      endelse
   endfor
endfor

map_e = double( fft( amn_e, /inverse, /double))
map_b = double( fft( amn_b, /inverse, /double))

end
