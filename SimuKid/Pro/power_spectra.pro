
whoami = 'nico' ; 'aina'

;; Outplot format
ps = 0
png= 0
nside = 1024 ; 512

if strupcase(whoami) eq 'AINA' then begin
   input_cl_file = "/Users/andriaai/NIKA/Processing/idl-libs/Healpix/test/cl.fits"
   map_dir = "/Users/andriaai/NIKA/Processing/Labtools/AA/Dev"
endif else begin
   input_cl_file = !nika.soft_dir+"/idl-libs/Healpix/test/cl.fits"
   map_dir = "/Data/Planck/Dust"

   ;; To have B modes in the CMB simulation
   ;;readcol, '/Data/CLS_totCls.dat',  l, cl1, cl2, cl3, cl4, format='DDDDD'
   readcol, "/Users/ponthien/Projects/NIKA/Soft/Papers_Proposal_Notes/SimuKid/Pro/cmb_totCls_r0.001.dat", $
            l, cl1, cl2, cl3, cl4, format='DDDDD'
   wind, 1, 1, /free
   plot, l, cl1, /xs, /ys, /xlog, yra=[1d-6,1d5], /ylog
   oplot, l, cl2, col=70
   oplot, l, cl3, col=250
   oplot, l, abs(cl4), col=150

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
endelse

npix = nside2npix( nside)

;; Theoretical Cl in microK^2
cl_in = mrdfits( input_cl_file, 1, cl_header)
ncl_in = n_elements(cl_in)

;; Produce maps
isynfast, input_cl_file, cmb_maps, nlmax=3*nside-1, nside=nside
;for i=0, 2 do mollview, cmb_maps[*,i]
;stop

;; Generate spurious signals 
map_x      = dblarr(npix,3)
map_x[*,0] = cmb_maps[*,0]^2
map_x[*,1] = 2.d0*cmb_maps[*,0]*cmb_maps[*,1]
map_x[*,2] = 2.d0*cmb_maps[*,0]*cmb_maps[*,2]

;; ;; Compute spurious power spectra
;; ispice, map_x, cover, lx, clxt, clxe, clxb, clxte
;; ;; remove monopole and dipole
;; lx    = lx[2:*]
;; clxt  = clxt[2:*]
;; clxe  = clxe[2:*]
;; clxb  = clxb[2:*]
;; clxte = clxte[2:*]

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

;; Pixel window functions
wlpix = mrdfits( !nika.soft_dir+"/idl-libs/Healpix/data/pixel_window_n"+$
                 string( nside, format='(I4.4)')+".fits", 1)
;; remove monopole and dipole
wlpix = wlpix[2:ncl_x-1]

;; Tolerance on BB contamination
tol = 1.d-3 ; microK^2 at l=100
tol = [1d-1,1d-2,1d-3]
ntol = n_elements(tol)

l_in      = dindgen(ncl_in)
fl_in     = l_in*(l_in+1)/(2*!dpi)

;; Cosmic Variance on BB
fsky = 1
c_var = sqrt(2./(fsky*(2*l_in+1.d0)))*cl_in.curl
tol = [1.d0, 0.1d0, 0.01d0]
ntol = n_elements(tol)

;;---------------------------------------------------------------
;; All Plots for quicklook
xra = [1, 5000]
cmb_col = !p.color
xcol = [70, 150, 250, 100, 200, 40]
wind, 1, 1, /free, /large
my_multiplot, 1, 2, pp, pp1, /rev
;; TT
yra = [10, max( reform( fl_x*clxt))]
plot, lx, fl_x*clxt, /xs, /ys, xra=xra, yra=yra, $
      xtitle='!12l', $
      ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n', $
      /ylog, /xlog, position=pp1[0,*], $
      title='Spurious Dust x Dust components'
oplot, l_in, fl_in*cl_in.temperature, col=0
oplot, lx, fl_x*clxt, col=xcol[0]
legendastro, ['CMB', 'X'], line=0, col=[cmb_col, xcol]

;; Polar
ymin = 1d-6
ymax = max([max(fl_x*clxe), max(fl_x*clxb), max( abs(fl_x*clxte))])
yra = [ymin, ymax]

wind, 1, 1, /free
plot, lx, fl_x*clxe, /xs, /ys, xra=xra, yra=yra, $
      xtitle='!12l', $
      ytitle='!12l(l+1)C!dl!n /2!7p!3 !7l!3K!u2!n', $
      /ylog, /xlog
oplot, l_in, fl_in*cl_in.gradient, col=cmb_col
oplot, l_in, fl_in*cl_in.curl, col=cmb_col
oplot, l_in, fl_in*abs(cl_in.G_T), col=cmb_col
oplot, lx, fl_x*clxe, col=xcol[1]
oplot, lx, fl_x*clxb, col=xcol[2]
oplot, lx, fl_x*abs(clxte), col=xcol[3]
for itol=0, ntol-1 do oplot, l_in, l_in*0 + fl_in*tol[itol]*c_var, line=2
legendastro, ['CMB', $
              '!7D!3EE', $
              '!7D!3BB', $
              '!7D!3TE'], line=0, $
             col=[cmb_col, xcol[1:3]], /bottom

;; Derive spec
w100 = where( l_in eq 100)
wx100 = where( lx eq 100)
epsilon = sqrt( tol*cl_in[w100].curl/(clxb[wx100])[0])
for itol=0, ntol-1 do print, ";; tol, epsilon: "+strtrim(tol[itol],2)+", "+strtrim(epsilon[itol],2)
;; tol, epsilon: 1.0000000, 0.00024614222
;; tol, epsilon: 0.10000000, 7.7837004e-05
;; tol, epsilon: 0.010000000, 2.4614222e-05

;; All spectra on the same plot
yra = [1d-8, 1d10]
wind, 1, 1, /free, /large
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
for itol=0, ntol-1 do oplot, lx, epsilon[itol]^2*fl_x*clxb, col=200
legendastro, ['CMB', $
              '!7D!3TT', $
              '!7D!3EE', $
              '!7D!3BB', $
              '!7D!3TE'], line=0, $
             col=[cmb_col, xcol[0:3]], /bottom


end
