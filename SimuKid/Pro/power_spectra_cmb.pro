
;; Outplot format
ps  = 0
png = 0

nside = 256 ; 1024 ; 512

map_dir = '$SK_DIR/Maps'

;; Input Cl's to have B modes in the CMB simulation
readcol, '$SK_DIR/Cl/cmb_totCls_r0.001.dat', $
         l, cl1, cl2, cl3, cl4, format='DDDDD'
wind, 1, 1, /free
plot,  l, cl1, /xs, /ys, /xlog, yra=[1d-6,1d5], /ylog
oplot, l, cl2, col=70
oplot, l, cl3, col=250
oplot, l, abs(cl4), col=150

;; Reshape for synfast:
;; Need to add monopole and dipole for synfast
cl_in = dblarr(max(l)+1,4)
;; prefactor included in the CAMB's output files
fl = l*(l+1)/(2*!dpi)
cl_in[2:*,0] = cl1/fl
cl_in[2:*,1] = cl2/fl
cl_in[2:*,2] = cl3/fl
cl_in[2:*,3] = cl4/fl

input_cl_file = 'mycl.fits'
cl2fits, cl_in, input_cl_file

;; Theoretical Cl in microK^2
cl_in = mrdfits( input_cl_file, 1, cl_header)
ncl_in = n_elements(cl_in)

;; Produce maps
isynfast, input_cl_file, cmb_maps, nlmax=3*nside-1, nside=nside
;for i=0, 2 do mollview, cmb_maps[*,i]
;stop

;; Convenient variables
l_in      = dindgen(ncl_in)
fl_in     = l_in*(l_in+1)/(2*!dpi)

;; Generate spurious signals 
npix = nside2npix( nside)
map_x      = dblarr(npix,3)
map_x[*,0] = cmb_maps[*,0]^2
map_x[*,1] = 2.d0*cmb_maps[*,0]*cmb_maps[*,1]
map_x[*,2] = 2.d0*cmb_maps[*,0]*cmb_maps[*,2]

;; Power spectra of spurious components
ianafast, map_x, cl_x, /double, nlmax=3*nside-1, /polar, /ring
lx = dindgen(n_elements(cl_x[*,0]))
lx = lx[2:*]
clxt  = reform( cl_x[2:*,0])
clxe  = reform( cl_x[2:*,1])
clxb  = reform( cl_x[2:*,2])
clxte = reform( cl_x[2:*,3])
clxtb = reform( cl_x[2:*,4])
clxeb = reform( cl_x[2:*,5])
ncl_x = n_elements(lx)
fl_x  = lx*(lx+1)/(2.*!dpi)

;; All spectra on the same plot
yra = [1d-7, 1d9]
xra = [2, 3000]
thick=2
xra = [1, 5000]
cmb_col = !p.color
xcol = [70, 150, 250, 100, 200, 40]

if ps eq 0 then wind, 1, 1, /free, /large
outplot, file='cmb_power_spectra', ps=ps, png=png, thick=thick
plot, lx, fl_x*clxt, /xs, /ys, xra=xra, yra=yra, $
      xtitle='!12l', $
      ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n', $
      /ylog, /xlog
oplot, l_in, fl_in*cl_in.temperature, col=0
oplot, lx, fl_x*clxt, col=xcol[0]
oplot, l_in, fl_in*cl_in.gradient, col=cmb_col
oplot, l_in, fl_in*cl_in.curl, col=cmb_col
oplot, l_in, fl_in*abs(cl_in.G_T), col=cmb_col
oplot, lx, fl_x*clxe, col=xcol[1]
oplot, lx, fl_x*clxb, col=xcol[2]
oplot, lx, fl_x*abs(clxte), col=xcol[3]
xyouts, 2.2, (fl_in*cl_in.temperature)[2]*1.3, 'TT'
xyouts, 2.2, (fl_in*cl_in.gradient)[2]*1.3, 'EE'
xyouts, 2.2, (fl_in*cl_in.curl)[2]*1.3, 'BB (r=0.001)'
xyouts, 2.2, (fl_in*abs(cl_in.g_t))[2]*1.3, 'TE'
legendastro, ['!7D!3TT', $
              '!7D!3EE', $
              '!7D!3BB', $
              '!7D!3TE'], line=0, $
             col=xcol[0:3]
xyouts, 10, 1d5, 'only cmb', orient=45, chars=2
outplot, /close, /verb

end
