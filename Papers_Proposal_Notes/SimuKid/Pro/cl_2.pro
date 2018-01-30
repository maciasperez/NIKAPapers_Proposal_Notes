;; from cl.pro
ps = 0 ;1
png= 0

;; Theory
input_cl_file = "/Users/andriaai/NIKA/Processing/idl-libs/Healpix/test/cl.fits"
cl_in = mrdfits( input_cl_file, 1)
ncl = n_elements(cl_in)

;; Generate a map
;; nside_th = 512
nside_th = 2048
isynfast, input_cl_file, map, nside=nside

;; Compute its power spectrum
ianafast, map, cl_out, /ring, /polar, nlmax=2*nside_th
nclout = n_elements(cl_out[*,0])
l_out = dindgen(nclout)
cl_out = cl_out[2:*,*]
l_out = l_out[2:*,*]

;; Pixel window functions
wlpix = mrdfits( "/Users/andriaai/NIKA/Processing/idl-libs/Healpix/data/pixel_window_n"+$
                 string(nside_th,format='(I4.4)')+".fits", 1)
wlpix = wlpix[2:nclout-1]

l = dindgen( ncl)
fl = l*(l+1)/(2*!dpi)
fl_out = l_out*(l_out+1)/(2.d0*!dpi)

;; Carte Planck
map_dir = "/Users/andriaai/NIKA/Processing/Labtools/AA/Dev"
map_planck = mrdfits( map_dir+"/HFI_SkyMap_353-field-IQU_2048_R2.02_full.fits", 1, header)
npix = n_elements(map_planck.i_stokes)
nside = npix2nside(npix)
ipnest = lindgen(npix)
nest2ring, nside, ipnest, ipring
;; bidon = dindgen(npix)
;; bidon[ipring] = map_planck.i_stokes
;; mollview, bidon 
map_i = dindgen(npix)
map_q = dindgen(npix)
map_u = dindgen(npix)
map_i[ipring] = map_planck.i_stokes*1000.d0
map_q[ipring] = map_planck.q_stokes*1000.d0
map_u[ipring] = map_planck.u_stokes*1000.d0

;; Tolerance on BB contamination
tol = 1.d-3 ; microK^2 at l=100
tol = [1d-1,1d-2,1d-3]
ntol = n_elements(tol)

;; ;Generate Cl with Planck map
map = dblarr(npix,3)
map[*,0] = map_i
map[*,1] = map_q
map[*,2] = map_u

ianafast, map, cl, /ring, /polar
cl = cl[2:*,*]

ncl_in    = n_elements(cl[*,0])
l_in      = dindgen(ncl_in)
fl_in     = l_in*(l_in+1)/(2*!dpi)
nclout_in = n_elements(cl)

;; Generate spurious signals 
map_x = dblarr(npix,3)
map_x[*,0] = map_i^2
map_x[*,1] = map_i * map_q
map_x[*,2] = map_i * map_u

ianafast, map_x, cl_planck_x, /ring, /polar;, nlmax=3*nside
cl_planck_x = cl_planck_x[2:*,*]

ncl    = n_elements(cl_planck_x[*,0])
l      = dindgen(ncl)
fl     = l*(l+1)/(2*!dpi)
nclout = n_elements(cl_planck_x)

w = where( l eq 100)
yra  = [1d-6,1d8]
ytol = [1d-3,1d-4,1d-5]
thick=2
epsilon_max = dblarr(ntol)
for i = 0,ntol-1 do begin
   epsilon_max[i] = (tol[i]*0.1/(fl_in*abs(cl_planck_x[*,3]))[w])[0]
endfor
print, 'r=', tol
print, 'epsilon^2 =', epsilon_max
print, 'epsilon =', sqrt(epsilon_max)
wind, 1, 1, /free, /large
outplot, file='cl2', png=png, ps=ps
plot, l, fl*cl_planck_x[*,0], /xs,/ys,xra=[1, max(l)], title='Spurious components', xtitle='!12l', $
      ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n',/ylog,/xlog, yra=yra,xthick=thick,ythick=thick,thick=thick
oplot, l, fl*cl_planck_x[*,1], col=250,thick=thick
oplot, l, fl*cl_planck_x[*,2], col=150,thick=thick
oplot, l, fl*abs(cl_planck_x[*,3]), col=70,thick=thick
for i = 0, ntol-1 do begin
   oplot, l, epsilon_max[i]*fl*cl_planck_x[*,2], col=150,thick=thick
   xyouts, 2, ytol[i],'r='+string(tol[i], format='(F5.3)')
endfor
legendastro, ['TT', 'EE','BB','TE'], line=0, col=[0,250,150,70], $
             box=0
outplot, /close, /verb


;; Meme figure mais sans TT et avec les val theoriques de Cl(BB), 24-10
yra  = [1d-10,1d5]
ytol = [1d-5,1d-6,1d-7]
wind, 1, 1, /free, /large
outplot, file='cl2_spurious', png=png, ps=ps
plot, l, fl*cl_planck_x[*,1], /xs,/ys,xra=[1, max(l)], title='Spurious components', xtitle='!12l', $
      ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n',/ylog,/xlog, yra=yra,xthick=thick,ythick=thick,thick=thick
oplot, l, fl*cl_planck_x[*,2], col=250,thick=thick
oplot, l, fl*abs(cl_planck_x[*,3]), col=70,thick=thick
for i = 0, ntol-1 do begin
   oplot, l_out, fl_out*epsilon_max[i]*abs(cl_in.g_t)/wlpix.temperature/wlpix.polarization, col=250,thick=thick
   xyouts, 2, ytol[i],'r='+string(tol[i], format='(F5.3)')
endfor
legendastro, ['EE','BB','TE'], line=0, col=[0,250,70], $
             box=0
outplot, /close, /verb


end
